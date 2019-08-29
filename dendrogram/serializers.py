from rest_framework import serializers
from dendrogram.models import Batch, Profile, Dendrogram


class BatchSerializer(serializers.ModelSerializer):
    class Meta:
        model = Batch
        fields = ('id', 'created', 'prof_num', 'linkage')


class ProfileSerializer(serializers.ModelSerializer):
    class Meta:
        model = Profile
        fields = ('id', 'batch_id', 'file')


class DendrogramSerializer(serializers.ModelSerializer):
    class Meta:
        model = Dendrogram
        fields = ('id', 'created', 'linkage', 'png_file', 'pdf_file', 'svg_file', 'newick_file', 'emf_file')
