docker run --name some-rabbit -d --env-file rabbitmq.env -p 5672:5672 -p 15672:15672 rabbitmq:3.7-management-alpine
docker run --name some-celery -d --env-file env.list -v /home/a504082002/Projects/Benga:/proj a504082002/benga:dev bash -c "cd /proj && python manage.py migrate && celery -A benga worker -l info"