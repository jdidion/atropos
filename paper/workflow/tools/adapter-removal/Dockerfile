# Here's how to build and deploy.
# Note: for me to be able to deploy this, it has to be tagged with my repo
# name. I'm not trying to take credit for anyone's work :)
# 1. docker build -f Dockerfile -t MikkelSchubert/adapterremoval:latest .
# 2. docker tag MikkelSchubert/adapterremoval:latest jdidion/adapterremoval
# 3. docker login -e <email> -u <user> -p <password>
# 4. docker push jdidion/adapterremoval

# Base Alpine
FROM alpine

# GCC
RUN apk add --no-cache g++ zlib-dev bzip2-dev

# Clone ngs-bits
RUN git clone --recursive https://github.com/MikkelSchubert/adapterremoval && \
    cd adapterremoval && \
    make 

ENTRYPOINT ["/adapterremoval/bin/AdapterRemoval"]
CMD ["--help"]
