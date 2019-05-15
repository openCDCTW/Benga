sudo apt-key adv --keyserver hkp://keyserver.ubuntu.com:80 --recv 9DA31620334BD75D9DCB49F368818C72E52529D4
echo "deb [ arch=amd64 ] https://repo.mongodb.org/apt/ubuntu bionic/mongodb-org/3.6 multiverse" \
    | sudo tee /etc/apt/sources.list.d/mongodb-org-3.6.list
sudo apt update
sudo apt install -y postgresql-client-10 mongodb-org-shell