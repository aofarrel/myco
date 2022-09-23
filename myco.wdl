version 1.0

import "https://raw.githubusercontent.com/aofarrel/clockwork-wdl/1.0.0/tasks/map_reads.wdl" as clckwrk_map_reads
import "https://raw.githubusercontent.com/aofarrel/clockwork-wdl/1.0.0/tasks/rm_contam.wdl" as clckwrk_rm_contam
import "https://raw.githubusercontent.com/aofarrel/clockwork-wdl/1.0.0/tasks/variant_call_one_sample.wdl" as clckwrk_var_call


task depth {
	input {
		File sam

		# runtime attributes
		Int addldisk = 250
		Int cpu      = 16
		Int retries  = 1
		Int memory   = 32
		Int preempt  = 1
	}
	String basestem_sam = basename(sam, ".sam")

	command <<<
	samtools sort -u ~{sam} > sorted_~{basestem_sam}.sam
	samtools depth -as ~{sam} > depth_~{basestem_sam}
	>>>

	runtime {
		cpu: cpu
		docker: "ashedpotatoes/iqbal-unofficial-clockwork-mirror:latest"
		disks: "local-disk " + finalDiskSize + " HDD"
		maxRetries: "${retries}"
		memory: "${memory} GB"
		preemptible: "${preempt}"
	}

	output {
		depth = glob("depth_*")[0]
	}

}


workflow myco {
	input:
		#File tarball_raw_ref
		#File tarball_decontaminated_ref
		#File tarball_H37Rv_ref
		#Array[String] samples
		#Array[Array[File]] fastqs
		Array[File] sam

	scatter(zip(samples, fastqs))
		call depth {
			input:
				sam = sam
		}
	
}

parameter_meta {
	"tarball_decontaminated_ref": "Indexed decontamination reference; output of clockwork reference_prepare"
	"tarball_H37Rv_ref": "Indexed H37Rv reference; output of clockwork reference_prepare"
	"fastqs": "An array of arrays. Each inner array represents one sample's FASTQs."
	"samples": "Sample names. Each sample corresponds to an inner array in the fastqs array at the same index."
}