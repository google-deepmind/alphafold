FROM nvidia/cuda:12.2.2-cudnn8-devel-ubuntu22.04
# FROM ubuntu:jammy

# ------------------------------------------------------------------------------
# Create a non-root user to use if preferred 
# - see https://aka.ms/vscode-remote/containers/non-root-user.
# ------------------------------------------------------------------------------
ARG USERNAME=vscode
ARG USER_UID=1000
ARG USER_GID=$USER_UID
# Create the user
RUN groupadd --gid $USER_GID $USERNAME \
    && useradd --uid $USER_UID --gid $USER_GID -m $USERNAME \
    # [Optional] Add sudo support. Omit if you don't need to install software after connecting.
    && apt-get update \
    && apt-get install -y sudo \
    && echo $USERNAME ALL=\(root\) NOPASSWD:ALL > /etc/sudoers.d/$USERNAME \
    && chmod 0440 /etc/sudoers.d/$USERNAME

# ------------------------------------------------------------------------------
# Install zsh and oh-my-zsh
# ------------------------------------------------------------------------------
RUN apt-get update && apt-get install -yq zsh sudo curl wget jq vim git-core gnupg locales && apt-get clean \
    && sh -c "$(curl -fsSL https://raw.githubusercontent.com/robbyrussell/oh-my-zsh/master/tools/install.sh)" || true \
    && sudo chsh -s /bin/zsh

# ------------------------------------------------------------------------------
# Install miniconda
# ------------------------------------------------------------------------------
# credits: @pangyuteng
# refer to: https://gist.github.com/pangyuteng/f5b00fe63ac31a27be00c56996197597
# Use the above args during building https://docs.docker.com/engine/reference/builder/#understand-how-arg-and-from-interact
# Choices
# PLATFORM: Linux, MacOSX, or Windows
# OS_TYPE : x86_64, arm64, ppc64le, s390x, or x86, armv7l
ARG CONDA_VER=latest
ARG PLATFORM=Linux
ARG OS_TYPE=x86_64
# Install miniconda to /miniconda
RUN curl -LO "http://repo.continuum.io/miniconda/Miniconda3-${CONDA_VER}-Linux-${OS_TYPE}.sh" \
    && bash Miniconda3-${CONDA_VER}-Linux-${OS_TYPE}.sh -p /miniconda -b \
    && rm Miniconda3-${CONDA_VER}-Linux-${OS_TYPE}.sh 
ENV PATH=/miniconda/bin:${PATH}
RUN conda update -y conda

# ------------------------------------------------------------------------------
# Install Docker, Docker-from-Docker 
# ------------------------------------------------------------------------------
# ref: https://github.com/microsoft/vscode-dev-containers/tree/main/containers/docker-from-docker 
# -------- Install Docker CE CLI --------
RUN apt-get update \
    && apt-get install -y apt-transport-https ca-certificates curl gnupg2 lsb-release \
    && curl -fsSL https://download.docker.com/linux/$(lsb_release -is | tr '[:upper:]' '[:lower:]')/gpg | apt-key add - 2>/dev/null \
    && echo "deb [arch=amd64] https://download.docker.com/linux/$(lsb_release -is | tr '[:upper:]' '[:lower:]') $(lsb_release -cs) stable" | tee /etc/apt/sources.list.d/docker.list \
    && apt-get update \
    && apt-get install -y docker-ce-cli

# -------- Install Docker Compose --------
RUN LATEST_COMPOSE_VERSION=$(curl -sSL "https://api.github.com/repos/docker/compose/releases/latest" | grep -o -P '(?<="tag_name": ").+(?=")') \
    && curl -sSL "https://github.com/docker/compose/releases/download/${LATEST_COMPOSE_VERSION}/docker-compose-$(uname -s)-$(uname -m)" -o /usr/local/bin/docker-compose \
    && chmod +x /usr/local/bin/docker-compose 
# NOTE: Remember add this following to your devcontainer.json: 
# "mounts": [ "source=/var/run/docker.sock,target=/var/run/docker.sock,type=bind" ]

# -------- Enable non-root user, vscode, to run docker --------
ARG NONROOT_USER=vscode

# Default to root only access to the Docker socket, set up non-root init script
RUN touch /var/run/docker-host.sock \
    && ln -s /var/run/docker-host.sock /var/run/docker.sock \
    && apt-get update \
    && apt-get -y install socat

# Create docker-init.sh to spin up socat
RUN echo "#!/bin/sh\n\
    sudoIf() { if [ \"\$(id -u)\" -ne 0 ]; then sudo \"\$@\"; else \"\$@\"; fi }\n\
    sudoIf rm -rf /var/run/docker.sock\n\
    ((sudoIf socat UNIX-LISTEN:/var/run/docker.sock,fork,mode=660,user=${NONROOT_USER} UNIX-CONNECT:/var/run/docker-host.sock) 2>&1 >> /tmp/vscr-docker-from-docker.log) & > /dev/null\n\
    \"\$@\"" >> /usr/local/share/docker-init.sh \
    && chmod +x /usr/local/share/docker-init.sh

# Cleanup
RUN apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* \
    && conda clean -a -y && pip cache purge

# ********************************************************
# * Anything else you want to do like clean up goes here *
# ********************************************************
# [Optional] Set the default user. Omit if you want to keep the default as root.
# If specified non-root user, add sudo when running e.g. apt-get
USER $USERNAME
RUN sudo usermod -aG sudo $USERNAME

# # create conda env
# COPY create-env.sh /home/vscode/create-env.sh
# RUN sudo chmod +x /home/vscode/create-env.sh && \
#     sudo chown vscode:vscode /home/vscode/create-env.sh && \
#     zsh /home/vscode/create-env.sh && \
#     rm /home/vscode/create-env.sh

# ------------------------------------------------------------------------------
# Install zsh and oh-my-zsh
# ------------------------------------------------------------------------------
RUN sh -c "$(curl -fsSL https://raw.githubusercontent.com/robbyrussell/oh-my-zsh/master/tools/install.sh)" || true \
    && sudo chsh -s /bin/zsh

# ------------------------------------------------------------------------------
# End of base image
# ------------------------------------------------------------------------------
WORKDIR /home/vscode

# VS Code overrides ENTRYPOINT and CMD when executing `docker run` by default.
# Setting the ENTRYPOINT to docker-init.sh will configure non-root access to
# the Docker socket if "overrideCommand": false is set in devcontainer.json.
# The script will also execute CMD if you need to alter startup behaviors.
ENTRYPOINT [ "/usr/local/share/docker-init.sh" ]
CMD [ "sleep", "infinity" ]