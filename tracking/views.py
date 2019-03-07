import json
from rest_framework import generics
from rest_framework.views import APIView
from rest_framework.response import Response
from rest_framework import status
from django.http import Http404
from tracking.serializers import ProfileSerializer, TrackedResultsSerializer,\
    TrackingSerializer
from tracking.tasks import track
from tracking.models import Profile, TrackedResults


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


class TrackedResultsList(generics.ListCreateAPIView):
    queryset = TrackedResults.objects.all()
    serializer_class = TrackedResultsSerializer


class TrackedResultsDetail(APIView):
    def get_object(self, pk):
        try:
            return TrackedResults.objects.get(pk=pk)
        except TrackedResults.DoesNotExist:
            raise Http404

    def get(self, request, pk, format=None):
        results = self.get_object(pk)
        serializer = TrackedResultsSerializer(results)
        data = serializer.data.copy()
        data["json"] = json.loads(results.json.read().decode("utf-8"))
        return Response(data)

    def put(self, request, pk, format=None):
        results = self.get_object(pk)
        serializer = TrackedResultsSerializer(results, data=request.data)
        if serializer.is_valid():
            serializer.save()
            return Response(serializer.data)
        return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)

    def delete(self, request, pk, format=None):
        results = self.get_object(pk)
        results.delete()
        return Response(status=status.HTTP_204_NO_CONTENT)


class Tracking(APIView):
    def post(self, request, format=None):
        serializer = TrackingSerializer(data=request.data)
        if serializer.is_valid():
            track.delay(str(serializer.data["id"]), str(serializer.data["profile_db"]))
            return Response(serializer.data, status=status.HTTP_202_ACCEPTED)
        return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)
