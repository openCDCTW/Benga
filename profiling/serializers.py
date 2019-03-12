from rest_framework import serializers
from profiling.models import Batch, Sequence, Profile


class BatchSerializer(serializers.ModelSerializer):
    class Meta:
        model = Batch
        fields = ('id', 'created')


class SequenceSerializer(serializers.ModelSerializer):
    class Meta:
        model = Sequence
        fields = ('id', 'batch_id', 'file', 'occurrence', 'database')


class ProfileSerializer(serializers.ModelSerializer):
    class Meta:
        model = Profile
        fields = ('id', 'created', 'zip', 'occurrence', 'database')
