from rest_framework import serializers
from profiling.models import UploadBatch, Sequence, Profile, Dendrogram


class UploadBatchSerializer(serializers.ModelSerializer):
    class Meta:
        model = UploadBatch
        fields = ('id', 'created')


class SequenceSerializer(serializers.ModelSerializer):
    class Meta:
        model = Sequence
        fields = ('id', 'batch_id', 'file')


class ProfileSerializer(serializers.ModelSerializer):
    class Meta:
        model = Profile
        fields = ('id', 'created', 'file', 'occurrence', 'database')


class DendrogramSerializer(serializers.ModelSerializer):
    class Meta:
        model = Dendrogram
        fields = ('id', 'created', 'png_file', 'pdf_file', 'svg_file', 'newick_file')
