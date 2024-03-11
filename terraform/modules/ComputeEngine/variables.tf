variable "project_id"{
    type = string
    description = "project ID"

}

variable "instance_name"{
    type = string
    description = "This is the name of the compute engine instance "
    default = "insmed-first-instance"
}

variable "zone"{
    type = string
    description = "This is the name of the zone "
    default = "us-central1-a"
}

variable "machine_type"{
    type = string
    description = "This is the type of machine  "
    default = "n2-standard-2"
}
variable "image"{
    type = string
    description = "Image for the Disk"
}

variable "network_name" {
  type = string 
  description = "name of the vpc in which vm is running"

}

variable "subnet_name" {
  type = string 
  description = "name of the subnetwork in which vm is running"
  
}

variable "labels" {
  type = map(string)
  description = "Labels asscociated for the instance created"
}


variable "startup_script" {
  type = string
  description = "Startup scripts which runs when vm starts "
}

# variable "guest_accelerator" {
#   type = object({
#     type = string
#     count = number
#   })
#   description = "hardware accelerators that can be attached to a virtual machine (VM) instance to enhance its performance for specific workloads"

# }

variable "accelerator_type" {
  description = "Type of guest accelerator"
  type        = string
}

variable "accelerator_count" {
  description = "Number of guest accelerators"
  type        = number
  default     = 1
}

variable "service_account_email" {
  description = "Email of the service account to be attached to the VM"
  type = string
}