# Low Resource Mode
Low Resource Mode runs the decontamination and variant calling tasks with barebones compute resources. This is designed for two use cases:
1) People who are running on a shared institute compute who don't want to be rude
2) People who are running on the cloud and want to gamble that this will be cheaper

The reason why (2) is a gamble for GCP (including Terra) users is that decontamination and variant calling both default to setting the WDL runtime attribute preemptible to 1. This means that the tasks will attempt to run in a preemptible instance once, and if the task fails, it will try again on a non-preemptible machine. Preemptible instances can save you >50% on cloud costs if they don't get preempted, so anything that meaningfully increases the chances they will be preempted is a bit of a gamble. If your samples are pretty small or you've enabled downsampling, I have a hunch it'll be worth it, but I have not tested this extensively and offer no guarantees.

AWS has something similar in concept to preemptible machines, although I don't know if they get called by the WDL runtime preemptible variable.