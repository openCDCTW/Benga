from rest_framework import serializers
from tracking.models import Sequence, TrackedResults


class SequenceSerializer(serializers.ModelSerializer):
    class Meta:
        model = Sequence
        fields = ('id', 'created', 'file')


class TrackedResultsSerializer(serializers.ModelSerializer):
    class Meta:
        model = TrackedResults
        fields = ('id', 'json')


class TrackingSerializer(serializers.ModelSerializer):
    class Meta:
        model = TrackedResults
        fields = ('id', 'database')
