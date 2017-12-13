FROM debian:latest

RUN apt update && apt upgrade -y && \
    apt install gfortran cmake python3 git -y && \
    git clone https://github.com/njoy/NJOY2016.git &&\
    cd NJOY2016 && mkdir build && cd build &&\
    cmake .. && \
    make

ENTRYPOINT [ "/bin/bash" ]
