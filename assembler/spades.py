#!/usr/bin/env python

import os
import shutil
import sys

spades_home = os.path.abspath(sys.path[0])

if os.path.realpath(__file__) == "/usr/bin/spades.py":
    spades_home = "/usr/share/spades/"

sys.path.append(os.path.join(spades_home, "src/tools/spades_pipeline/"))
sys.path.append(os.path.join(spades_home, "src/tools/quality/"))
sys.path.append(os.path.join(spades_home, "src/tools/quality/libs"))

import support
from process_cfg import *

def error(err_str, prefix="== Error == "):
    print("\n\n" + prefix + " " + err_str + "\n\n")
    exit(1)

def warning(warn_str, prefix="== Warning == "):
    print("\n\n" + prefix + " " + warn_str + "\n\n")

def prepare_config_bh(filename, cfg):

    subst_dict = dict()
    cfg.working_dir = os.path.abspath(cfg.working_dir)

    import glob
    input_reads = []
    if not type(cfg.input_reads) is list:
        cfg.input_reads = [cfg.input_reads]
    for read in cfg.input_reads:
        input_reads.extend(glob.glob(os.path.abspath(os.path.expandvars(read))))

    if len(input_reads) == 0:
        error("The given wildcards do not correspond to existing files: " + str(cfg.input_reads))

    cfg.input_reads = input_reads

    if len(cfg.input_reads) == 1:
        import bh_aux
        cfg.input_reads = bh_aux.split_paired_file(cfg)

    subst_dict["input_numfiles"]           = len(cfg.input_reads)

    for i, input_read in enumerate(cfg.input_reads):
        subst_dict["input_file_" + str(i)]  = input_read

    subst_dict["input_gzipped"]             = bool_to_str(cfg.input_reads[0].endswith(".gz"))
    subst_dict["input_working_dir"]         = cfg.working_dir
    subst_dict["general_max_iterations"]    = cfg.max_iterations
    subst_dict["general_max_nthreads"]      = cfg.max_threads
    subst_dict["count_merge_nthreads"]      = cfg.max_threads
    subst_dict["bayes_nthreads"]            = cfg.max_threads
    subst_dict["expand_nthreads"]           = cfg.max_threads
    subst_dict["correct_nthreads"]          = cfg.max_threads
    subst_dict["general_hard_memory_limit"] = cfg.max_memory

    if cfg.__dict__.has_key("qvoffset"):
        subst_dict["input_qvoffset"]        = cfg.qvoffset

    substitute_params(filename, subst_dict) 

def prepare_config_spades(filename, cfg, prev_K, last_one):

    subst_dict = dict()
    cfg.working_dir = os.path.abspath(cfg.working_dir)

    subst_dict["dataset"]            = cfg.dataset
    subst_dict["output_base"]        = cfg.working_dir
    subst_dict["additional_contigs"] = cfg.additional_contigs
    subst_dict["entry_point"]        = 'construction'
    subst_dict["developer_mode"]     = str(cfg.developer_mode).lower()
    subst_dict["SAM_writer_enable"]  = bool_to_str(cfg.generate_sam_files)

    if last_one:
        subst_dict["gap_closer_enable"] = "true"
        if cfg.paired_mode:
            subst_dict["paired_mode"] = "true"
    else:
        subst_dict["paired_mode"] = "false"
        subst_dict["gap_closer_enable"] = "false"

    if prev_K:
        subst_dict["use_additional_contigs"] = "false"
    else:
        subst_dict["use_additional_contigs"] = "true"

    substitute_params(filename, subst_dict)

def check_config(cfg, config_filename):

    if (not cfg.has_key("error_correction")) and (not cfg.has_key("assembly")):
        error("wrong config! You should specify either 'error_correction' section (for reads error correction) or 'assembly' one (for assembling) or both!")
        return False

    if not cfg["common"].__dict__.has_key("output_dir"):
        error("wrong config! You should specify output_dir!")
        return False

    if not cfg["common"].__dict__.has_key("output_to_console"):
        cfg["common"].__dict__["output_to_console"] = True

    if not cfg["common"].__dict__.has_key("developer_mode"):
        cfg["common"].__dict__["developer_mode"] = False

    if not cfg["common"].__dict__.has_key("project_name"):
        cfg["common"].__dict__["project_name"] = os.path.splitext(os.path.basename(config_filename))[0]

    if not cfg["common"].__dict__.has_key("compilation_dir"):
        cfg["common"].__dict__["compilation_dir"] = os.path.join(os.getenv('HOME'), '.spades/precompiled/')
    else:
        cfg["common"].compilation_dir = os.path.abspath(cfg["common"].compilation_dir)

    cfg["common"].output_dir = os.path.join(os.path.abspath(os.path.expandvars(cfg["common"].output_dir)), cfg["common"].project_name)

    return True

def main():

    CONFIG_FILE = ""
    if os.path.isfile("spades_config.info"):
        CONFIG_FILE = "spades_config.info"
    elif os.path.isfile(os.path.join(spades_home, "spades_config.info")):
        CONFIG_FILE = os.path.join(spades_home, "spades_config.info")
    if len(sys.argv) > 1 and os.path.isfile(sys.argv[1]):
       CONFIG_FILE = sys.argv[1]
    if not CONFIG_FILE:
        print >> sys.stderr, "Usage : python ", sys.argv[0], " <config file>"
        return

    print("Using config file: " + CONFIG_FILE)
    os.environ["cfg"] = os.path.dirname(os.path.abspath(CONFIG_FILE))

    cfg = load_config_from_info_file(CONFIG_FILE)

    if not check_config(cfg, CONFIG_FILE):
        return

    bh_dataset_filename = ""
    if cfg.has_key("error_correction"):
        bh_cfg = cfg["error_correction"]

        if not bh_cfg.__dict__.has_key("output_dir"):
            bh_cfg.__dict__["output_dir"] = os.path.join(os.path.expandvars(cfg["common"].output_dir), "corrected")
        else:
            bh_cfg.__dict__["output_dir"] = os.path.expandvars(bh_cfg.output_dir)
            warning("output_dir (" + bh_cfg.output_dir + ") will be used for error correction instead of the common one (" + cfg["common"].output_dir + ")")
        if not bh_cfg.__dict__.has_key("dataset_name"):
            bh_cfg.__dict__["dataset_name"] = "dataset"
        
        bh_cfg = merge_configs(cfg["error_correction"], cfg["common"])
        
        bh_cfg.__dict__["working_dir"] = os.path.join(bh_cfg.output_dir, "tmp")

        start_bh = True
        if os.path.exists(bh_cfg.output_dir):
            bh_dataset_filename = os.path.abspath(os.path.join(bh_cfg.output_dir, bh_cfg.dataset_name + ".info"))            
            if os.path.exists(bh_dataset_filename):
                question = ["WARNING! It looks like error correction was already done!", 
                            "Dataset " + bh_dataset_filename + " already exists!",
                            "Do you want to overwrite this dataset and start error correction again?"]
                answer = support.question_with_timer(question, 10, 'n')
                if answer == 'n':
                    start_bh = False
                    print("\n===== Error correction skipped\n")
                else:
                    os.remove(bh_dataset_filename)
            else:
                bh_dataset_filename = ""
            #    shutil.rmtree(bh_cfg.output_dir)

        if start_bh:
            if not os.path.exists(bh_cfg.working_dir):
                os.makedirs(bh_cfg.working_dir)

            log_filename = os.path.join(bh_cfg.output_dir, "correction.log")
            bh_cfg.__dict__["log_filename"] = log_filename

            shutil.copy(CONFIG_FILE, bh_cfg.output_dir)

            print("\n===== Error correction started. Log can be found here: " + bh_cfg.log_filename + "\n")

            tee = support.Tee(log_filename, 'w', console=bh_cfg.output_to_console)

            err_code = 0
            try:
                bh_dataset_filename = run_bh(bh_cfg)
            except support.spades_error as err:
                print err.err_str
                err_code = err.code

            tee.free()

            print("\n===== Error correction finished. Log can be found here: " + bh_cfg.log_filename + "\n")
            if err_code:
               exit(err_code)

    if cfg.has_key("assembly"):
        spades_cfg = merge_configs(cfg["assembly"], cfg["common"])        
        if not spades_cfg.__dict__.has_key("generate_sam_files"):
            spades_cfg.__dict__["generate_sam_files"] = False

        if bh_dataset_filename != "":
            if spades_cfg.__dict__.has_key("dataset"):
                warning("dataset created during error correction (" + bh_dataset_filename + ") will be used instead of the one from config file (" + spades_cfg.dataset + ")!")
            spades_cfg.__dict__["dataset"] = bh_dataset_filename
        spades_cfg.dataset = os.path.abspath(os.path.expandvars(spades_cfg.dataset))
        if not os.path.isfile(spades_cfg.dataset):
            error("dataset " + spades_cfg.dataset + " doesn't exist!")
            exit(1)
        if not check_dataset(spades_cfg.dataset):
            error("" "incorrect dataset " + spades_cfg.dataset + ": list of files with reads should be arounded with double quotes!")
            exit(1)

        def make_working_dir(output_dir):
            import datetime
            name = "spades_" + datetime.datetime.now().strftime("%m.%d_%H.%M.%S")
            working_dir = os.path.join(output_dir, name)
            os.makedirs(working_dir)
            latest = os.path.join(output_dir, "latest")
            if os.path.islink(latest):
                os.remove(latest)
            if not os.path.exists(latest):
                os.symlink(name, latest)
            return working_dir

        spades_cfg.__dict__["working_dir"] = make_working_dir(spades_cfg.output_dir)
        spades_cfg.__dict__["log_filename"] = os.path.join(spades_cfg.working_dir, "assembly.log")
        spades_cfg.__dict__["result_contigs"] = os.path.join(spades_cfg.working_dir, spades_cfg.project_name + ".fasta")
        spades_cfg.__dict__["additional_contigs"] = os.path.join(spades_cfg.working_dir, "simplified_contigs.fasta")

        shutil.copy(CONFIG_FILE, spades_cfg.working_dir)
        shutil.copy(spades_cfg.dataset, spades_cfg.working_dir)

        print("\n===== Assembling started. Log can be found here: " + spades_cfg.log_filename + "\n")

        tee = support.Tee(spades_cfg.log_filename, 'w', console=spades_cfg.output_to_console)
        
        err_code = 0
        try:
            if cfg.has_key("quality_assessment"):
                run_spades(spades_cfg, cfg["quality_assessment"])
            else:
                run_spades(spades_cfg)
        except support.spades_error as err:
            print err.err_str
            err_code = err.code

        tee.free()

        print("\n===== Assembling finished. Log can be found here: " + spades_cfg.log_filename + "\n")
        if err_code:
            exit(err_code)

    print("\n===== SPAdes pipeline finished\n") 

def run_bh(cfg):

    dst_configs = os.path.join(cfg.output_dir, "configs")
    if os.path.exists(dst_configs):
        shutil.rmtree(dst_configs)
    shutil.copytree(os.path.join(spades_home, "configs"), dst_configs)
    cfg_file_name = os.path.join(dst_configs, "hammer", "config.info")

    prepare_config_bh(cfg_file_name, cfg)

    import build
    build.build_hammer(cfg, spades_home)

    execution_home = os.path.join(cfg.compilation_dir, 'build_hammer')
    command = os.path.join(execution_home, "hammer", "hammer") + " " + os.path.abspath(cfg_file_name)

    print("\n== Running error correction tool: " + command + "\n")
    support.sys_call(command)

    import bh_aux
    dataset_str = bh_aux.generate_dataset(cfg, cfg.working_dir, cfg.input_reads)    
    dataset_filename = os.path.abspath(os.path.join(cfg.output_dir, cfg.dataset_name + ".info"))
    dataset_file = open(dataset_filename, "w")
    dataset_file.write(dataset_str)
    dataset_file.close()
    print("\n== Dataset created: " + dataset_filename + "\n")

    shutil.rmtree(cfg.working_dir)

    return dataset_filename

def run_spades(cfg, quality_cfg = None):

    if type(cfg.iterative_K) is int:
        cfg.iterative_K = [cfg.iterative_K]
    cfg.iterative_K = sorted(cfg.iterative_K)

    import build
    build.build_spades_n_copy(cfg, spades_home)

    count = 0
    prev_K = None

    for K in cfg.iterative_K:
        count += 1

        dst_configs = os.path.join(cfg.working_dir, "config_K" + str(K), "configs")
        shutil.copytree(os.path.join(spades_home, "configs"), dst_configs)
        cfg_file_name = os.path.join(dst_configs, "debruijn", "config.info")

        prepare_config_spades(cfg_file_name, cfg, prev_K, count == len(cfg.iterative_K))
        prev_K = K

        execution_home = os.path.join(cfg.compilation_dir, 'build' + str(K))
        command = os.path.join(execution_home, "debruijn", "spades") + " " + os.path.abspath(cfg_file_name)

        print("\n== Running assembler: " + command + "\n")

        support.sys_call(command, execution_home)

        #dataset_id = os.path.splitext(os.path.basename(cfg.dataset))[0]
        latest = os.path.join(cfg.working_dir, "K%d" % (K), "latest")
        latest = os.readlink(latest)
        latest = os.path.join(cfg.working_dir, "K%d" % (K), latest)
        os.symlink(os.path.relpath(latest, cfg.working_dir), os.path.join(cfg.working_dir, "link_K%d" % (K)))

    shutil.copyfile(os.path.join(latest, "final_contigs.fasta"), cfg.result_contigs)
    os.remove(cfg.additional_contigs)

    if quality_cfg:
        print("\n== Running quality assessment tools: " + cfg.log_filename + "\n")

        args = [cfg.result_contigs]
        #dataset_filename = cfg.dataset
        #dataset = load_config_from_info_file(dataset_filename)["common"]

        if quality_cfg.__dict__.has_key("reference"):
            args.append("-R")
            args.append(os.path.abspath(os.path.expandvars(quality_cfg.reference)) )
        if quality_cfg.__dict__.has_key("genes"):
            args.append("-G")
            args.append(os.path.abspath(os.path.expandvars(quality_cfg.genes)) )
        if quality_cfg.__dict__.has_key("operons"):
            args.append("-O")
            args.append(os.path.abspath(os.path.expandvars(quality_cfg.operons)) )
        quality_output_dir = os.path.join(cfg.working_dir, "quality_results")
        args.append("-o")
        args.append(quality_output_dir)
        import quality
        quality.main(args, lib_dir=os.path.join(spades_home, "src/tools/quality/libs"))

    print ""
    print "All the resulting information can be found here: " + cfg.working_dir
    print " * Resulting contigs are called " + os.path.basename(cfg.result_contigs)
    if quality_cfg:
        print " * Assessment of their quality is in " + quality_output_dir + "/"
    print ""
    print "Thank you for using SPAdes!"

if __name__ == '__main__':
    main()
