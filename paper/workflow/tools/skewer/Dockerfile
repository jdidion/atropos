# Here's how to build and deploy.
# Note: for me to be able to deploy this, it has to be tagged with my repo
# name. I'm not trying to take credit for anyone's work :)
# 1. docker build -f Dockerfile -t replimoc/skewer:latest .
# 2. docker tag replimoc/skewer:latest jdidion/skewer
# 3. docker login -e <email> -u <user> -p <password>
# 4. docker push jdidion/skewer

# Base Alpine
FROM alpine

# GCC
RUN apk add --no-cache g++

# Clone ngs-bits
RUN git clone --recursive https://github.com/relipmoc/skewer && \
    cd skewer && \
    make 

ENTRYPOINT ["/skewer/bin/skewer"]
CMD ["--help"]
