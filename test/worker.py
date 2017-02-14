import pika
import time

connection = pika.BlockingConnection(pika.ConnectionParameters(host='rabbit'))
channel = connection.channel()

queue = 'task_queue'

channel.queue_declare(queue=queue, durable=True)
print(' [*] Waiting for messages. To exit press CTRL+C')


def callback(ch, method, properties, body):
    print(" [x] Received %r" % (body,))
    time.sleep(body.count(b'.'))
    print(" [x] Done")
    ch.basic_ack(delivery_tag=method.delivery_tag)


channel.basic_qos(prefetch_count=1)
channel.basic_consume(callback, queue=queue)

channel.start_consuming()

