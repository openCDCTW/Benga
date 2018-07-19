from rest_framework import serializers


class EchoSerializer(serializers.Serializer):
    number = serializers.IntegerField()
