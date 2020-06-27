from rest_framework import serializers
from core.models import Arquivo

class CoreSerializers(serializers.ModelSerializer):
    class Meta:
        model = Arquivo
        depth = 1
        fields = ['id','files']