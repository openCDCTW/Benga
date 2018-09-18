from django.urls import path

from . import consumers

websocket_urlpatterns = [
    path('ws/profiling/<uuid:pk>/$', consumers.ProfilingConsumer),
]