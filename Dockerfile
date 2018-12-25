FROM a504082002/seqtool
MAINTAINER Yueh-Hua Tu <a504082002@gmail.com>

# Install dependencies
RUN pip install --trusted-host pypi.python.org -r requirements.txt --upgrade

# Copy project
WORKDIR /benga
COPY . /benga

# Initialize database tables
RUN python3 manage.py makemigrations && \
    python3 manage.py migrate

CMD [".", "run.sh"]
