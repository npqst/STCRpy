# Code adapted from https://github.com/lucidrains/egnn-pytorch/blob/main/egnn_pytorch/egnn_pytorch_geometric.py


# MIT License

# Copyright (c) 2021 Phil Wang, Eric Alcaide

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.


import torch
from torch import nn, einsum, broadcast_tensors
import torch.nn.functional as F

from einops import rearrange, repeat
from einops.layers.torch import Rearrange

# helper functions


def exists(val):
    return val is not None


def safe_div(num, den, eps=1e-8):
    res = num.div(den.clamp(min=eps))
    res.masked_fill_(den == 0, 0.0)
    return res


def batched_index_select(values, indices, dim=1):
    value_dims = values.shape[(dim + 1) :]
    values_shape, indices_shape = map(lambda t: list(t.shape), (values, indices))
    indices = indices[(..., *((None,) * len(value_dims)))]
    indices = indices.expand(*((-1,) * len(indices_shape)), *value_dims)
    value_expand_len = len(indices_shape) - (dim + 1)
    values = values[(*((slice(None),) * dim), *((None,) * value_expand_len), ...)]

    value_expand_shape = [-1] * len(values.shape)
    expand_slice = slice(dim, (dim + value_expand_len))
    value_expand_shape[expand_slice] = indices.shape[expand_slice]
    values = values.expand(*value_expand_shape)

    dim += value_expand_len
    return values.gather(dim, indices)


def fourier_encode_dist(x, num_encodings=4, include_self=True):
    x = x.unsqueeze(-1)
    device, dtype, orig_x = x.device, x.dtype, x
    scales = 2 ** torch.arange(num_encodings, device=device, dtype=dtype)
    x = x / scales
    x = torch.cat([x.sin(), x.cos()], dim=-1)
    x = torch.cat((x, orig_x), dim=-1) if include_self else x
    return x


def embedd_token(x, dims, layers):
    stop_concat = -len(dims)
    to_embedd = x[:, stop_concat:].long()
    for i, emb_layer in enumerate(layers):
        # the portion corresponding to `to_embedd` part gets dropped
        x = torch.cat([x[:, :stop_concat], emb_layer(to_embedd[:, i])], dim=-1)
        stop_concat = x.shape[-1]
    return x


# swish activation fallback


class Swish_(nn.Module):
    def forward(self, x):
        return x * x.sigmoid()


SiLU = nn.SiLU if hasattr(nn, "SiLU") else Swish_

# helper classes

# this follows the same strategy for normalization as done in SE3 Transformers
# https://github.com/lucidrains/se3-transformer-pytorch/blob/main/se3_transformer_pytorch/se3_transformer_pytorch.py#L95


class CoorsNorm(nn.Module):
    def __init__(self, eps=1e-8, scale_init=1.0):
        super().__init__()
        self.eps = eps
        scale = torch.zeros(1).fill_(scale_init)
        self.scale = nn.Parameter(scale)

    def forward(self, coors):
        norm = coors.norm(dim=-1, keepdim=True)
        normed_coors = coors / norm.clamp(min=self.eps)
        return normed_coors * self.scale


# global linear attention


class Attention(nn.Module):
    def __init__(self, dim, heads=8, dim_head=64):
        super().__init__()
        inner_dim = heads * dim_head
        self.heads = heads
        self.scale = dim_head**-0.5

        self.to_q = nn.Linear(dim, inner_dim, bias=False)
        self.to_kv = nn.Linear(dim, inner_dim * 2, bias=False)
        self.to_out = nn.Linear(inner_dim, dim)

    def forward(self, x, context, mask=None):
        h = self.heads

        q = self.to_q(x)
        kv = self.to_kv(context).chunk(2, dim=-1)

        q, k, v = map(lambda t: rearrange(t, "b n (h d) -> b h n d", h=h), (q, *kv))
        dots = einsum("b h i d, b h j d -> b h i j", q, k) * self.scale

        if exists(mask):
            mask_value = -torch.finfo(dots.dtype).max
            mask = rearrange(mask, "b n -> b () () n")
            dots.masked_fill_(~mask, mask_value)

        attn = dots.softmax(dim=-1)
        out = einsum("b h i j, b h j d -> b h i d", attn, v)

        out = rearrange(out, "b h n d -> b n (h d)", h=h)
        return self.to_out(out)


class GlobalLinearAttention(nn.Module):
    def __init__(self, *, dim, heads=8, dim_head=64):
        super().__init__()
        self.norm_seq = nn.LayerNorm(dim)
        self.norm_queries = nn.LayerNorm(dim)
        self.attn1 = Attention(dim, heads, dim_head)
        self.attn2 = Attention(dim, heads, dim_head)

        self.ff = nn.Sequential(
            nn.LayerNorm(dim),
            nn.Linear(dim, dim * 4),
            nn.GELU(),
            nn.Linear(dim * 4, dim),
        )

    def forward(self, x, queries, mask=None):
        res_x, res_queries = x, queries
        x, queries = self.norm_seq(x), self.norm_queries(queries)

        induced = self.attn1(queries, x, mask=mask)
        out = self.attn2(x, induced)

        x = out + res_x
        queries = induced + res_queries

        x = self.ff(x) + x
        return x, queries


# classes


class EGNN(nn.Module):
    def __init__(
        self,
        dim,
        edge_dim=0,
        m_dim=16,
        fourier_features=0,
        num_nearest_neighbors=0,
        dropout=0.0,
        init_eps=1e-3,
        norm_feats=False,
        norm_coors=False,
        norm_coors_scale_init=1e-2,
        update_feats=True,
        update_coors=True,
        only_sparse_neighbors=False,
        valid_radius=float("inf"),
        m_pool_method="sum",
        soft_edges=False,
        coor_weights_clamp_value=None,
        return_edges=False,
    ):
        super().__init__()
        assert m_pool_method in {
            "sum",
            "mean",
        }, "pool method must be either sum or mean"
        assert (
            update_feats or update_coors
        ), "you must update either features, coordinates, or both"

        self.fourier_features = fourier_features

        edge_input_dim = (fourier_features * 2) + (dim * 2) + edge_dim + 1
        dropout = nn.Dropout(dropout) if dropout > 0 else nn.Identity()

        self.edge_mlp = nn.Sequential(
            nn.Linear(edge_input_dim, edge_input_dim * 2),
            dropout,
            SiLU(),
            nn.Linear(edge_input_dim * 2, m_dim),
            SiLU(),
        )

        self.edge_gate = (
            nn.Sequential(nn.Linear(m_dim, 1), nn.Sigmoid()) if soft_edges else None
        )

        self.node_norm = nn.LayerNorm(dim) if norm_feats else nn.Identity()
        self.coors_norm = (
            CoorsNorm(scale_init=norm_coors_scale_init) if norm_coors else nn.Identity()
        )

        self.m_pool_method = m_pool_method

        self.node_mlp = (
            nn.Sequential(
                nn.Linear(dim + m_dim, dim * 2),
                dropout,
                SiLU(),
                nn.Linear(dim * 2, dim),
            )
            if update_feats
            else None
        )

        self.coors_mlp = (
            nn.Sequential(
                nn.Linear(m_dim, m_dim * 4), dropout, SiLU(), nn.Linear(m_dim * 4, 1)
            )
            if update_coors
            else None
        )

        self.num_nearest_neighbors = num_nearest_neighbors
        self.only_sparse_neighbors = only_sparse_neighbors
        self.valid_radius = valid_radius

        self.coor_weights_clamp_value = coor_weights_clamp_value

        self.init_eps = init_eps
        self.apply(self.init_)

        self.return_edges = return_edges

    def init_(self, module):
        if type(module) in {nn.Linear}:
            # seems to be needed to keep the network from exploding to NaN with greater depths
            nn.init.normal_(module.weight, std=self.init_eps)

    def forward(self, feats, coors, edges=None, mask=None, adj_mat=None):
        (
            b,
            n,
            d,
            device,
            fourier_features,
            num_nearest,
            valid_radius,
            only_sparse_neighbors,
        ) = (
            *feats.shape,
            feats.device,
            self.fourier_features,
            self.num_nearest_neighbors,
            self.valid_radius,
            self.only_sparse_neighbors,
        )

        if exists(mask):
            num_nodes = mask.sum(dim=-1)

        use_nearest = num_nearest > 0 or only_sparse_neighbors

        rel_coors = rearrange(coors, "b i d -> b i () d") - rearrange(
            coors, "b j d -> b () j d"
        )
        rel_dist = (rel_coors**2).sum(dim=-1, keepdim=True)

        i = j = n

        if use_nearest:
            ranking = rel_dist[..., 0].clone()

            if exists(mask):
                rank_mask = mask[:, :, None] * mask[:, None, :]
                ranking.masked_fill_(~rank_mask, 1e5)

            if exists(adj_mat):
                if len(adj_mat.shape) == 2:
                    adj_mat = repeat(adj_mat.clone(), "i j -> b i j", b=b)

                if only_sparse_neighbors:
                    num_nearest = int(adj_mat.float().sum(dim=-1).max().item())
                    valid_radius = 0

                self_mask = rearrange(
                    torch.eye(n, device=device, dtype=torch.bool), "i j -> () i j"
                )

                adj_mat = adj_mat.masked_fill(self_mask, False)
                ranking.masked_fill_(self_mask, -1.0)
                ranking.masked_fill_(adj_mat, 0.0)

            nbhd_ranking, nbhd_indices = ranking.topk(
                num_nearest, dim=-1, largest=False
            )

            nbhd_mask = nbhd_ranking <= valid_radius

            rel_coors = batched_index_select(rel_coors, nbhd_indices, dim=2)
            rel_dist = batched_index_select(rel_dist, nbhd_indices, dim=2)

            if exists(edges):
                edges = batched_index_select(edges, nbhd_indices, dim=2)

            j = num_nearest

        if fourier_features > 0:
            rel_dist = fourier_encode_dist(rel_dist, num_encodings=fourier_features)
            rel_dist = rearrange(rel_dist, "b i j () d -> b i j d")

        if use_nearest:
            feats_j = batched_index_select(feats, nbhd_indices, dim=1)
        else:
            feats_j = rearrange(feats, "b j d -> b () j d")

        feats_i = rearrange(feats, "b i d -> b i () d")
        feats_i, feats_j = broadcast_tensors(feats_i, feats_j)

        edge_input = torch.cat((feats_i, feats_j, rel_dist), dim=-1)

        if exists(edges):
            edge_input = torch.cat((edge_input, edges), dim=-1)

        m_ij = self.edge_mlp(edge_input)

        if exists(self.edge_gate):
            m_ij = m_ij * self.edge_gate(m_ij)

        if exists(mask):
            mask_i = rearrange(mask, "b i -> b i ()")

            if use_nearest:
                mask_j = batched_index_select(mask, nbhd_indices, dim=1)
                mask = (mask_i * mask_j) & nbhd_mask
            else:
                mask_j = rearrange(mask, "b j -> b () j")
                mask = mask_i * mask_j

        if exists(self.coors_mlp):
            coor_weights = self.coors_mlp(m_ij)
            coor_weights = rearrange(coor_weights, "b i j () -> b i j")

            rel_coors = self.coors_norm(rel_coors)

            if exists(mask):
                coor_weights.masked_fill_(~mask, 0.0)

            if exists(self.coor_weights_clamp_value):
                clamp_value = self.coor_weights_clamp_value
                coor_weights.clamp_(min=-clamp_value, max=clamp_value)

            coors_out = (
                einsum("b i j, b i j c -> b i c", coor_weights, rel_coors) + coors
            )
        else:
            coors_out = coors

        if exists(self.node_mlp):
            if exists(mask):
                m_ij_mask = rearrange(mask, "... -> ... ()")
                m_ij = m_ij.masked_fill(~m_ij_mask, 0.0)

            if self.m_pool_method == "mean":
                if exists(mask):
                    # masked mean
                    mask_sum = m_ij_mask.sum(dim=-2)
                    m_i = safe_div(m_ij.sum(dim=-2), mask_sum)
                else:
                    m_i = m_ij.mean(dim=-2)

            elif self.m_pool_method == "sum":
                m_i = m_ij.sum(dim=-2)

            normed_feats = self.node_norm(feats)
            node_mlp_input = torch.cat((normed_feats, m_i), dim=-1)
            node_out = self.node_mlp(node_mlp_input) + feats
        else:
            node_out = feats

        if self.return_edges:
            if exists(num_nearest):
                num_neighbours = num_nearest
            else:
                num_neighbours = n
            edges_out = torch.zeros(
                (b, n, n, m_ij.shape[-1])
            )  # initialise full edge matrix
            edges_out.scatter_(
                2, nbhd_indices.unsqueeze(-1).expand(-1, -1, -1, m_ij.shape[-1]), m_ij
            )
            # assert torch.stack([edges_out[0, i, nbhd_indices[0, i, idx]] == m_ij[0, i, idx] for idx in range(num_neighbours) for i in range(n)]).all()          # check assignemnt has worked, comment out at runtime
            return node_out, coors_out, edges_out

        return node_out, coors_out
