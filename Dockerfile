FROM a504082002/biopython-mkl
MAINTAINER Yueh-Hua Tu <a504082002@gmail.com>

# Copy project
WORKDIR /benga
COPY . /benga

# Install backend dependencies
RUN pip install --trusted-host pypi.python.org -r requirements.txt --upgrade

# Install frontend dependencies
RUN curl -o- https://raw.githubusercontent.com/creationix/nvm/v0.33.11/install.sh | bash && \
    nvm install v10.13.0 && \
    npm install

CMD ["/bin/sh"]
