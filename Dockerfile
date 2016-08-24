FROM a504082002/seqtool-python3
MAINTAINER a504082002 <a504082002@gmail.com>

RUN mkdir /program
ADD src /program

CMD ["/bin/bash"]

