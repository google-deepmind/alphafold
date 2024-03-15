
project_id = "quantiphi-sandox-env"
ComputeEngine = [{
  instance_name = "VM_NAME" # this will be replaced in by the github actions workflow
  machine_type  = "custom-8-53248" #Custom machine types can be formatted as custom-NUMBER_OF_CPUS-AMOUNT_OF_MEMORY_MB 52GB => 53248MB
  zone          = "us-central1-a"
  image = "projects/ml-images/global/images/c0-deeplearning-common-cu113-v20230925-debian-10"
  labels = {
    environment = "dev"
  }
  network_name      = "vpc-githubactions-test"
  subnet_name       = "sbn-githubactions-test-us-central1"
  accelerator_type  = "nvidia-tesla-v100"
  accelerator_count = 1
  service_account_email = "557417765420-compute@developer.gserviceaccount.com"
  startup_script = <<-EOF
        #!/bin/bash
          sudo /opt/deeplearning/install-driver.sh
          EOF
}]


