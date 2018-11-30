from django.http import Http404
from rest_framework import generics, status
from rest_framework.response import Response
from rest_framework.views import APIView

from dendrogram.models import Profile, Dendrogram
from dendrogram.serializers import ProfileSerializer, DendrogramSerializer, PlotingSerializer
from dendrogram.tasks import plot_dendrogram


class ProfileList(generics.ListCreateAPIView):
    queryset = Profile.objects.all()
    serializer_class = ProfileSerializer


class ProfileDetail(APIView):
    def get_object(self, pk):
        try:
            return Profile.objects.get(pk=pk)
        except Profile.DoesNotExist:
            raise Http404

    def get(self, request, pk, format=None):
        sequence = self.get_object(pk)
        serializer = ProfileSerializer(sequence)
        return Response(serializer.data)

    def put(self, request, pk, format=None):
        sequence = self.get_object(pk)
        serializer = ProfileSerializer(sequence, data=request.data)
        if serializer.is_valid():
            serializer.save()
            return Response(serializer.data)
        return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)

    def delete(self, request, pk, format=None):
        sequence = self.get_object(pk)
        sequence.delete()
        return Response(status=status.HTTP_204_NO_CONTENT)


class DendrogramList(generics.ListCreateAPIView):
    queryset = Dendrogram.objects.all()
    serializer_class = DendrogramSerializer


class DendrogramDetail(APIView):
    def get_object(self, pk):
        try:
            return Dendrogram.objects.get(pk=pk)
        except Profile.DoesNotExist:
            raise Http404

    def get(self, request, pk, format=None):
        sequence = self.get_object(pk)
        serializer = DendrogramSerializer(sequence)
        return Response(serializer.data)

    def put(self, request, pk, format=None):
        sequence = self.get_object(pk)
        serializer = DendrogramSerializer(sequence, data=request.data)
        if serializer.is_valid():
            serializer.save()
            return Response(serializer.data)
        return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)

    def delete(self, request, pk, format=None):
        sequence = self.get_object(pk)
        sequence.delete()
        return Response(status=status.HTTP_204_NO_CONTENT)


class Plotting(APIView):
    def post(self, request, format=None):
        serializer = PlotingSerializer(data=request.data)
        if serializer.is_valid():
            plot_dendrogram.s(str(serializer.data["id"]))
            return Response(serializer.data, status=status.HTTP_202_ACCEPTED)
        return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)
