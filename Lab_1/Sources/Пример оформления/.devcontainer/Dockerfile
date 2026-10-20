FROM ghcr.io/astral-sh/uv:0.10.0 AS uv_binary
FROM mcr.microsoft.com/devcontainers/python:1-3.13-bookworm
COPY --from=uv_binary /uv /uvx /usr/local/bin/
