FROM a504082002/biopython-mkl
MAINTAINER Yueh-Hua Tu <a504082002@gmail.com>

# Copy project
WORKDIR /benga
COPY . /benga

# Install dependencies
RUN pip install --trusted-host pypi.python.org -r requirements.txt --upgrade

CMD ["/bin/sh"]
