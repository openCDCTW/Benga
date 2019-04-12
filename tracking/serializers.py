from rest_framework import serializers
from tracking.models import Profile, TrackedResults


class ProfileSerializer(serializers.ModelSerializer):
    class Meta:
        model = Profile
        fields = ('id', 'created', 'file', 'profile_db')


class TrackedResultsSerializer(serializers.ModelSerializer):
    class Meta:
        model = TrackedResults
        fields = ('id', 'json')
