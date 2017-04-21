#################################################################
# Dockerfile
#
# Description:      Base image for data containers.
# Provides:         wget
# Base Image:       alpine
# Build Cmd:        docker build -f Dockerfile -t jdidion/data_base:latest .
#################################################################
FROM alpine
RUN apk --update add ca-certificates openssl wget && \
    update-ca-certificates && \
    rm -rf /var/cache/apk/*