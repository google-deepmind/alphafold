import requests
from google.cloud import compute_v1

METADATA_URL = "http://metadata.google.internal/computeMetadata/v1/"
METADATA_HEADERS = {"Metadata-Flavor": "Google"}

project_id = requests.get(url=METADATA_URL+"project/project-id",headers=METADATA_HEADERS).text
print("project_id",project_id)

zone = requests.get(url=METADATA_URL+"instance/zone",headers=METADATA_HEADERS).text
print("zone",zone)

instance_name = requests.get(url=METADATA_URL+"instance/name",headers=METADATA_HEADERS).text
print("instance name", instance_name)

compute_client = compute_v1.InstancesClient()
operation = compute_client.delete(project=project_id,zone=zone.split("/")[-1],instance=instance_name)
operation.result()
print("VM is deleted")
