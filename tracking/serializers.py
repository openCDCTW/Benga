from rest_framework import serializers
from tracking.models import Profile, TrackedResults, Tracking


class ProfileSerializer(serializers.ModelSerializer):
    class Meta:
        model = Profile
        fields = ('id', 'created', 'file')


class TrackedResultsSerializer(serializers.ModelSerializer):
    class Meta:
        model = TrackedResults
        fields = ('id', 'json')


class TrackingSerializer(serializers.ModelSerializer):
    class Meta:
        model = Tracking
        fields = ('id', 'allele_db', 'profile_db')
