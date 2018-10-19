from django.urls import path
from profiling.consumers import ProfilingConsumer

websocket_urlpatterns = [
    path('ws/profiling/<uuid:batch_id>/', ProfilingConsumer),
]