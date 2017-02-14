import pika
import sys

connection = pika.BlockingConnection(pika.ConnectionParameters(host='rabbit'))
channel = connection.channel()

queue = 'task_queue'

channel.queue_declare(queue=queue, durable=True)

message = ' '.join(sys.argv[1:]) or "Hello World!"
channel.basic_publish(exchange='', routing_key='hello', body=message)
print(" [x] Sent %r" % (message,))

connection.close()

