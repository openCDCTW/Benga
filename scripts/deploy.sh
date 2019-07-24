# backend
sudo apt update
sudo apt install -y git python3-pip virtualenv postgresql-client-10 rabbitmq-server
git clone https://github.com/openCDCTW/Benga.git
cd Benga/
git checkout develop
virtualenv venv -p python3
. venv/bin/activate
pip install -r requirements.txt --upgrade

# build frontend
curl -o- https://raw.githubusercontent.com/creationix/nvm/v0.33.11/install.sh | bash
export NVM_DIR="/root/.nvm"
[ -s "$NVM_DIR/nvm.sh" ] && \. "$NVM_DIR/nvm.sh"
[ -s "$NVM_DIR/bash_completion" ] && \. "$NVM_DIR/bash_completion"
# disconnect and reconnect
cd Benga/
. venv/bin/activate

##
nvm install v10.13.0
npm install
export NODE_OPTIONS=--max_old_space_size=4096
npm run build
# modify environment variables
. scripts/run_envs.sh  envs/api-server.env
. scripts/run_envs.sh  envs/nosql.env
. scripts/run_envs.sh  envs/sql.env
. scripts/generate_secret_key.sh

python manage.py check
python manage.py makemigrations
python manage.py migrate

sudo mkdir -p /var/www/benga/cgMLST/
sudo chown ubuntu /var/www/benga/cgMLST/
python manage.py collectstatic

sudo chown -R www-data /var/www/benga
sudo chmod -R 755 /var/www/benga

# configure
sudo systemctl restart benga.service