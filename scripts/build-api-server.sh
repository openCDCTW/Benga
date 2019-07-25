# Install dependencies
sudo apt-key adv --keyserver hkp://keyserver.ubuntu.com:80 --recv 9DA31620334BD75D9DCB49F368818C72E52529D4
echo "deb [ arch=amd64 ] https://repo.mongodb.org/apt/ubuntu bionic/mongodb-org/4.0 multiverse" \
    | sudo tee /etc/apt/sources.list.d/mongodb-org-4.0.list
sudo apt update
sudo apt install -y git python3 virtualenv postgresql-client-10 mongodb-org-shell mongo-tools rabbitmq-server libreoffice ncbi-blast+

# Deploy Benga project
git clone https://github.com/openCDCTW/Benga.git
cd Benga/
virtualenv venv -p python3
. venv/bin/activate
pip install -r requirements.txt
# modify env variables

# Restore postgresql
createdb -h postgresql.cc22wzgngpjq.us-east-1.rds.amazonaws.com --port=5432 -U centrallab --password ${DBNAME}
xz -d dump_Vibrio_cholerae_15-05-2019_10_21_11.xz | psql -h postgresql.cc22wzgngpjq.us-east-1.rds.amazonaws.com --port=5432 -U centrallab --password ${DBNAME}

# Restore mongodb
mongorestore --host 172.31.2.104:27017 -u centrallab -p ${PASSWORD} dump_${DBNAME}_*

# Setup rabbitmq
sudo rabbitmq-plugins enable rabbitmq_management
sudo rabbitmqctl add_user ${RABBITMQ_DEFAULT_USER} ${RABBITMQ_DEFAULT_PASS}
sudo rabbitmqctl add_vhost benga
sudo rabbitmqctl set_permissions -p benga ${RABBITMQ_DEFAULT_USER} ".*" ".*" ".*"