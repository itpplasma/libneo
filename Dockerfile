FROM debian:bookworm-slim

ADD setup/debian.sh /tmp/setup.sh
RUN sh /tmp/setup.sh
