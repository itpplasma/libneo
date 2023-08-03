FROM debian:bookworm-slim

ADD install_deps.sh /
RUN sh /install_deps.sh
