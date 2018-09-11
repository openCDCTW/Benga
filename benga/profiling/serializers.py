from rest_framework import serializers
from profiling.models import UploadBatch, Sequence, Profile


class UploadBatchSerializer(serializers.ModelSerializer):
    class Meta:
        model = UploadBatch
        fields = ('id', 'created')


class SequenceSerializer(serializers.ModelSerializer):
    class Meta:
        model = Sequence
        fields = ('id', 'batch_id', 'filename', 'file')


class ProfileSerializer(serializers.ModelSerializer):
    class Meta:
        model = Profile
        fields = ('id', 'created', 'file', 'occurrence', 'database')

    def create(self, validated_data):
        validated_data["file"] = validated_data["file"].read()
        print(validated_data)
        return Profile.objects.create(**validated_data)