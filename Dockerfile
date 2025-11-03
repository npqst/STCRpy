FROM python:3.12-slim

# STCRpy Docker Image with full analysis and visualization support
# Includes: STCRpy, PLIP (interaction profiling), PyMOL (3D visualization)
# Install system dependencies including OpenBabel and PyMOL dependencies
# Using python3-openbabel from Debian avoids building from source
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y \
    build-essential \
    wget \
    git \
    libxml2-dev \
    libxslt1-dev \
    zlib1g-dev \
    openbabel \
    libopenbabel7 \
    libopenbabel-dev \
    python3-openbabel \
    # PyMOL build dependencies
    libglew-dev \
    libpng-dev \
    libfreetype6-dev \
    libmsgpack-dev \
    python3-dev \
    libglm-dev \
    # Qt5 dependencies for PyMOL GUI
    libqt5core5a \
    libqt5gui5 \
    libqt5widgets5 \
    libqt5opengl5 \
    qt5-qmake \
    qtbase5-dev \
    libxcb-xinerama0 \
    libxkbcommon-x11-0 \
    && rm -rf /var/lib/apt/lists/*

# Install uv
COPY --from=ghcr.io/astral-sh/uv:latest /uv /usr/local/bin/uv

# Set working directory
WORKDIR /app

# Copy project files
COPY . .

# Create virtual environment with system-site-packages to access python3-openbabel
RUN uv venv /opt/venv --system-site-packages && \
    . /opt/venv/bin/activate && \
    uv pip install -e . && \
    uv pip install pymol-open-source PyQt5 && \
    ANARCI --build_models

# Install PLIP source code directly (avoids pip build issues)
# PLIP is pure Python and uses system openbabel bindings
RUN git clone https://github.com/pharmai/plip.git /opt/plip && \
    ln -s /opt/plip/plip /opt/venv/lib/python3.12/site-packages/plip

# Set environment to use venv
ENV PATH="/opt/venv/bin:$PATH"
ENV VIRTUAL_ENV="/opt/venv"

# Default command
CMD ["/bin/bash"]
