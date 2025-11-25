FROM python:3.12-slim

# STCRpy Docker Image with full analysis, ML, and visualization support
# Includes: STCRpy, PLIP, PyMOL, scikit-learn, PyTorch, transformers
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

WORKDIR /app
COPY . .

# Install STCRpy and dependencies
RUN pip install --no-cache-dir --root-user-action ignore -e ".[ml_datasets]" \
    && pip install --no-cache-dir --root-user-action ignore einops pymol-open-source PyQt5 \
    && ANARCI --build_models

# Install PLIP source code directly (avoids pip build issues)
# PLIP is pure Python and uses system openbabel bindings
RUN git clone https://github.com/pharmai/plip.git /opt/plip && \
    ln -s /opt/plip/plip /usr/local/lib/python3.12/site-packages/

CMD ["/bin/bash"]
