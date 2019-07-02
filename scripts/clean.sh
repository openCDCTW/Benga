sudo systemctl stop benga.service
sudo systemctl stop celery.service

sudo > /var/log/celery/worker.log
sudo > /var/log/uwsgi/benga.log

sudo rm /var/www/benga/cgMLST/media/*

dropdb -h postgresql.cc22wzgngpjq.us-east-1.rds.amazonaws.com --port=5432 -U centrallab --password benga
createdb -h postgresql.cc22wzgngpjq.us-east-1.rds.amazonaws.com --port=5432 -U centrallab --password benga

rm dendrogram/migrations/0001_initial.py
rm profiling/migrations/0001_initial.py
rm tracking/migrations/0001_initial.py

python manage.py makemigrations
python manage.py migrate

sudo systemctl start celery.service
sudo systemctl start benga.service