FROM a504082002/biopython-mkl
MAINTAINER Yueh-Hua Tu <a504082002@gmail.com>

# Copy project
WORKDIR /benga
COPY . /benga

# Install backend dependencies
RUN pip install --trusted-host pypi.python.org -r requirements.txt --upgrade

# Install frontend dependencies
RUN curl -o- https://raw.githubusercontent.com/creationix/nvm/v0.33.11/install.sh | bash && \
    export NVM_DIR="/root/.nvm" && \
    [ -s "$NVM_DIR/nvm.sh" ] && \. "$NVM_DIR/nvm.sh" && \
    [ -s "$NVM_DIR/bash_completion" ] && \. "$NVM_DIR/bash_completion" && \
    nvm install v10.13.0 && \
    npm install && \
    export NODE_OPTIONS=--max_old_space_size=4096 && \
    npm run build

CMD ["gunicorn", "benga.wsgi:application", "-w", "2", "-k", "uvicorn.workers.UvicornWorker"]
