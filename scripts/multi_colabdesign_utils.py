#============================================================================
# Multi ColabDesign Utils
#============================================================================
# Inlined verbatim from functions/multi_colabdesign_utils.py; imports hoisted to the cell above.

def binder_hallucination_multi(template_information, design_paths, design_name, advanced_settings, length, seed, chain, helicity_value, design_models):
    model_pdb_path = os.path.join(design_paths["Trajectory"], design_name+".pdb")
    # clear GPU memory for new trajectory
    clear_mem()
    passing = True

    print("updates")
    af_model_optimization = mk_afdesign_model(
        protocol="binder", 
        debug=False, 
        data_dir=advanced_settings["af_params_dir"], 
        use_multimer=advanced_settings["use_multimer_design"], 
        num_recycles=advanced_settings["num_recycles_design"],
        best_metric='loss'
    )

    # check for hotspots and initialize the model with corresponding hotsxpots:
    for template in template_information:
        if template_information[template]["template_hostspot_residues"] == "":
            template_information[template]["template_hostspot_residues"] = None
        af_model_optimization.prep_inputs(
                pdb_filename=template_information[template]["template_pdb"], 
                chain=chain, 
                binder_len=length, 
                hotspot=template_information[template]["template_hostspot_residues"], 
                seed=seed, 
                rm_aa=advanced_settings["omit_AAs"],
                # rm_template_ic=True
                )

        template_information[template]["batch"] = copy.deepcopy(af_model_optimization._inputs['batch'])
        template_information[template]["hotspot"] = copy.deepcopy(af_model_optimization.opt.get("hotspot"))
    
    ### Update weights based on specified settings
    af_model_optimization.opt["weights"].update({
        "pae":advanced_settings["weights_pae_intra"],
        "plddt":advanced_settings["weights_plddt"],
        "i_pae":advanced_settings["weights_pae_inter"],
        "con":advanced_settings["weights_con_intra"],
        "i_con":advanced_settings["weights_con_inter"],
    })

    # redefine intramolecular contacts (con) and intermolecular contacts (i_con) definitions
    af_model_optimization.opt["con"].update({
        "num":advanced_settings["intra_contact_number"],
        "cutoff":advanced_settings["intra_contact_distance"],
        "binary":False,
        "seqsep":9})
    af_model_optimization.opt["i_con"].update({
        "num":advanced_settings["inter_contact_number"],
        "cutoff":advanced_settings["inter_contact_distance"],
        "binary":False})

    ### additional loss functions: 
    # radius of gyration loss
    if advanced_settings["use_rg_loss"]:
        add_rg_loss(af_model_optimization, advanced_settings["weights_rg"])
    # interface pTM loss
    if advanced_settings["use_i_ptm_loss"]:
        add_i_ptm_loss(af_model_optimization, advanced_settings["weights_iptm"])
    # add target plddt loss 
    if advanced_settings["use_target_plddt_loss"]:
        add_target_plddt_loss(af_model_optimization, advanced_settings["weights_target_plddt"])
    # add the helicity loss
    add_helix_loss(af_model_optimization, helicity_value)


    # calculate the number of mutations to do based on the length of the protein
    greedy_tries = math.ceil(length * (advanced_settings["greedy_percentage"] / 100))

#SECTIO II-III run trajectory optimization
    ### start design algorithm based on selection
    if advanced_settings["design_algorithm"] == '2stage':
        # uses gradient descend to get a PSSM profile and then uses PSSM to bias the sampling of random mutations to decrease loss
        raise NotImplementedError("not yet implemented")

    elif advanced_settings["design_algorithm"] == '3stage':
        # 3 stage design using logits, softmax, and one hot encoding
        raise NotImplementedError("not yet implemented")

    elif advanced_settings["design_algorithm"] == 'greedy':
        # design by using random mutations that decrease loss
        raise NotImplementedError("not yet implemented")

    elif advanced_settings["design_algorithm"] == 'mcmc':
        # design by using random mutations that decrease loss
        raise NotImplementedError("not yet implemented")
    
    elif advanced_settings["design_algorithm"] == '4stage_duo':
        if len(template_information) == 0:
            print("Preparing Logits")
            raise NotImplementedError() 

        else:
            if advanced_settings["buffer_iterations"] > 0:
                print("Preparing Logits on one target trajectory")
                first_target_template_information = {k: v for k, v in [next(((k, v) for k, v in template_information.items() if v['template_targeting_strategy'] in ['', 'target']), (None, None))] if k}
                design_multi(af_model_optimization, first_target_template_information, advanced_settings, iters=advanced_settings["buffer_iterations"], e_soft=0, models=design_models, num_models=1, sample_models=advanced_settings["sample_models"], save_best=True, dropout=False, binder_length=length)
            
            if advanced_settings["use_loss_decay"]:
                print(f"Stage 1 adv.: Increasing the detarget loss over the next {advanced_settings['decay_iters']} iters")
                design_multi(af_model_optimization, template_information, advanced_settings, iters=advanced_settings["decay_iters"], e_soft=0, models=design_models, num_models=1, sample_models=advanced_settings["sample_models"], save_best=True, dropout=False, binder_length=length, decay_iters=advanced_settings["decay_iters"])
                
                print(f"Stage 1: Test Logits")
                design_multi(af_model_optimization, template_information, advanced_settings, iters=advanced_settings["test_iterations"], e_soft=0.9, models=design_models, num_models=1, sample_models=advanced_settings["sample_models"], save_best=True, dropout=False, binder_length=length, decay_scale=advanced_settings["default_detargeting_decay"])

            else:
                print("Stage 1: Test Logits")
                # optional raw-logits warm-up (e_soft=0, so soft stays 0 -- no softmax) before the soft ramp
                _warmup = int(advanced_settings.get("multi_stage1_warmup_iters", 0) or 0)
                if _warmup > 0:
                    design_multi(af_model_optimization, template_information, advanced_settings, iters=_warmup, e_soft=0, models=design_models, num_models=1, sample_models=advanced_settings["sample_models"], save_best=True, dropout=False, binder_length=length, decay_scale=advanced_settings["default_detargeting_decay"])
                design_multi(af_model_optimization, template_information, advanced_settings, iters=advanced_settings["test_iterations"], e_soft=0.9, models=design_models, num_models=1, sample_models=advanced_settings["sample_models"], save_best=True, dropout=False, binder_length=length, decay_scale=advanced_settings["default_detargeting_decay"])
            # determine pLDDT of best iteration according to lowest 'loss' value
            initial_plddt = get_best_plddt(af_model_optimization, length)

            # if best iteration has high enough confidence then continue
            if initial_plddt > 0.65:
                print("Initial trajectory pLDDT good, continuing: "+str(initial_plddt))
                if advanced_settings["optimise_beta"]:
                    # temporarily dump model to assess secondary structure
                    af_model_optimization.save_pdb(model_pdb_path)
                    _, beta, *_ = calc_ss_percentage(model_pdb_path, advanced_settings, 'B')
                    os.remove(model_pdb_path)

                    # if beta sheeted trajectory is detected then choose to optimise
                    if float(beta) > 15:
                        advanced_settings["soft_iterations"] = advanced_settings["soft_iterations"] + advanced_settings["optimise_beta_extra_soft"]
                        advanced_settings["temporary_iterations"] = advanced_settings["temporary_iterations"] + advanced_settings["optimise_beta_extra_temp"]
                        af_model_optimization.set_opt(num_recycles=advanced_settings["optimise_beta_recycles_design"])
                        print("Beta sheeted trajectory detected, optimising settings")

                # how many logit iterations left
                logits_iter = advanced_settings["soft_iterations"] - 50
                if logits_iter > 0:
                    print("Stage 1: Additional Logits Optimisation")
                    af_model_optimization.clear_best()
                    design_multi(af_model_optimization, template_information, advanced_settings, iters=logits_iter, e_soft=1, models=design_models, num_models=1, sample_models=advanced_settings["sample_models"], save_best=True, 
                                        dropout=False, ramp_recycles=False, binder_length=length)
                    af_model_optimization._tmp["seq_logits"] = af_model_optimization.aux["seq"]["logits"]
                    logit_plddt = get_best_plddt(af_model_optimization, length)
                    print("Optimised logit trajectory pLDDT: "+str(logit_plddt))
                else:
                    logit_plddt = initial_plddt

                # perform softmax trajectory design
                if advanced_settings["temporary_iterations"] > 0:
                    print("Stage 2: Softmax Optimisation")
                    af_model_optimization.clear_best()
                    design_multi(af_model_optimization, template_information, advanced_settings, advanced_settings["temporary_iterations"], e_temp=1e-2, models=design_models, num_models=1,
                                        sample_models=advanced_settings["sample_models"], ramp_recycles=False, save_best=True, dropout=False, soft=1, binder_length=length, decay_scale=advanced_settings["default_detargeting_decay"])
                    softmax_plddt = get_best_plddt(af_model_optimization, length)
                else:
                    softmax_plddt = logit_plddt

                # perform one hot encoding
                if softmax_plddt > 0.65:
                    print("Softmax trajectory pLDDT good, continuing: "+str(softmax_plddt))
                    if advanced_settings["hard_iterations"] > 0:
                        af_model_optimization.clear_best()
                        print("Stage 3: One-hot Optimisation")
                        design_multi(af_model_optimization, template_information, advanced_settings, advanced_settings["hard_iterations"], temp=1e-2, models=design_models, num_models=1,
                                        sample_models=advanced_settings["sample_models"], dropout=False, ramp_recycles=False, save_best=True,hard=1, soft=1, binder_length=length, decay_scale=advanced_settings["default_detargeting_decay"])
                        onehot_plddt = get_best_plddt(af_model_optimization, length)

                    if onehot_plddt > 0.65:
                        # perform greedy mutation optimisation
                        print("One-hot trajectory pLDDT good, continuing: "+str(onehot_plddt))
                        if advanced_settings["greedy_iterations"] > 0:
                            print("Stage 4: PSSM Semigreedy Optimisation")
                            af_model_optimization.clear_best()
                            design_semigreedy_multi(af_model_optimization, template_information, advanced_settings, iters=advanced_settings["greedy_iterations"], tries=greedy_tries, models=design_models, num_models=1, sample_models=advanced_settings["sample_models"], save_best=True, dropout=False, hard=1, soft=0, binder_length=length)
                            onehot_plddt = get_best_plddt(af_model_optimization, length)

                        if advanced_settings["mcmc_iterations"] > 0:
                            print("Stage 4: MCMC  Optimisation")
                            # design by using random mutations that decrease loss
                            half_life = round(advanced_settings["mcmc_iterations"] / 5, 0)
                            t_mcmc = 0.01
                            design_mcmc_multi(af_model_optimization, template_information, advanced_settings, steps=advanced_settings["mcmc_iterations"], half_life=half_life, T_init=t_mcmc, mutation_rate=1, seq_logits=None, sample_models=(advanced_settings["sample_models"]), save_best=True, verbose=True)


                    else:   
                        for template in template_information:
                            update_failures(template_information[template]["failure_csv"], 'Trajectory_one-hot_pLDDT')
                        af_model_optimization.aux['log']['terminate'] = "Skip"
                        print("One-hot trajectory pLDDT too low to continue: "+str(onehot_plddt))
                        passing=False
                        return passing, af_model_optimization

                else:
                    for template in template_information:
                        update_failures(template_information[template]["failure_csv"], 'Trajectory_softmax_pLDDT')
                    af_model_optimization.aux['log']['terminate'] = "Skip"
                    print("Softmax trajectory pLDDT too low to continue: "+str(softmax_plddt))
                    passing=False
                    return passing, af_model_optimization
            else:
                for template in template_information:
                    update_failures(template_information[template]["failure_csv"], 'Trajectory_logits_pLDDT')
                print("Initial trajectory pLDDT too low to continue: "+str(initial_plddt))
                af_model_optimization.aux['log']['terminate'] = "Skip"
                passing=False
                return passing, af_model_optimization
    else:
        print("ERROR: No valid design model selected")
        exit()

    af_model_optimization.set_seq()
    af_model_optimization.aux["log"]["terminate"] = ""
    binder_seq=af_model_optimization.get_seq(get_best=True)[0]
    trajectory_dir=design_paths["Trajectory"]

    for template in template_information:
        template_information[template]["trajectory_model"], template_information[template]["trajectory_pdb"] = process_trajectory_complex(
            binder_sequence=binder_seq, 
            design_name=design_name, 
            target_pdb=template_information[template]["template_pdb"], 
            chain=chain, 
            length=length, 
            trajectory_dir=trajectory_dir, 
            prediction_models=[0], 
            appendix=template_information[template]["template_name"], 
            failure_csv=template_information[template]["failure_csv"], 
            advanced_settings=advanced_settings
        )

    for template in template_information:
        if template_information[template]["trajectory_model"].aux["log"]["terminate"]:
            passing = False
            if template_information[template]["trajectory_model"].aux["log"]["terminate"] != "":
                shutil.move(template_information[template]["trajectory_pdb"], design_paths[f"Trajectory/{template_information[template]['trajectory_model'].aux['log']['terminate']}"])

    ### get the sampled sequence for plotting
    af_model_optimization.get_seqs()
    if advanced_settings["save_design_trajectory_plots"]:
        plot_trajectory(af_model_optimization, design_name, design_paths)

    return passing, template_information
        


def step_multi(model_to_optimize, template_information, advanced_settings, lr_scale=1.0, num_recycles=None, num_models=None, sample_models=None, models=None, backprop=True, callback=None, save_best=False, verbose=1, binder_length=None, iter_idx=None, decay_iters=None, decay_scale = None):
    # Multi-template aggregation controls, implemented as in
    # multi_colabdesign_utils_patched.py, including its defaults. Note those defaults
    # turn shuffling and RMS normalisation ON when the key is absent, so a settings
    # file that predates these keys changes behaviour; both settings files shipped
    # with this notebook set all five explicitly. Only step_multi (stages 1-3) honours
    # them -- the MCMC and semigreedy stages still aggregate the original way.
    shuffle_templates = advanced_settings.get("multi_shuffle_templates", True)
    grad_rms_norm     = advanced_settings.get("multi_grad_rms_norm", True)
    grad_rms_eps      = advanced_settings.get("multi_grad_rms_eps", 1e-8)
    dyn_beta          = advanced_settings.get("multi_dynamic_weight_beta", 0.0)
    bottleneck_alpha  = advanced_settings.get("multi_bottleneck_alpha", 0.0)

    template_order = list(template_information.keys())
    if shuffle_templates and len(template_order) > 1:
        template_order = [str(t) for t in np.random.permutation(template_order)]

    grads, base_ws = [], []
    # run optmiization for each template 
    for template in template_order:
        # empty & set template losses:
        template_information[template]["logging"] = {}

        # swap batch for current template:
        if template_information[template]["batch"] != None:
            model_to_optimize._inputs['batch'] = template_information[template]["batch"]
        # swap the hotspot too: prep_inputs() left opt["hotspot"] on whichever template was prepped last
        if template_information[template].get("hotspot") is None:
            model_to_optimize.opt.pop("hotspot", None)
        else:
            model_to_optimize.opt["hotspot"] = template_information[template]["hotspot"]

        # swap weights for targeting/detargeting:
        if template_information[template]["template_targeting_strategy"] == "target":
            model_to_optimize.opt["weights"].update(
                {"pae":(advanced_settings["weights_pae_intra"]   ),
                "plddt":(advanced_settings["weights_plddt"]      ),
                "i_pae":(advanced_settings["weights_pae_inter"]  ),
                "con":(advanced_settings["weights_con_intra"]    ),
                "i_con":(advanced_settings["weights_con_inter"]  ),
                "i_ptm":(advanced_settings["weights_iptm"]       ),
    })
        # scaling the weights according to detargeting scale
        elif template_information[template]["template_targeting_strategy"] == "detarget":
            model_to_optimize.opt["weights"].update(
                {"pae":(advanced_settings["weights_pae_intra"]   *(advanced_settings["detargeting_scaling_weights_pae_intra"])),
                "plddt":(advanced_settings["weights_plddt"]      *(advanced_settings["detargeting_scaling_weights_plddt"])),
                "i_pae":(advanced_settings["weights_pae_inter"]  *(advanced_settings["detargeting_scaling_weights_pae_inter"])),
                "con":(advanced_settings["weights_con_intra"]    *(advanced_settings["detargeting_scaling_weights_con_intra"])),
                "i_con":(advanced_settings["weights_con_inter"]  *(advanced_settings["detargeting_scaling_weights_con_inter"])),
                "i_ptm":(advanced_settings["weights_iptm"]       *(advanced_settings["detargeting_scaling_weights_iptm"])),
    })
            # Apply decay scaling if enabled
            if decay_iters is not None and iter_idx is not None:
                decay_scale = min(1.0, 0.1 + 0.9 * (iter_idx / decay_iters))
            # print("detargeting weights: ", model_to_optimize.opt["weights"])
            
        else:
            raise ValueError("No senseful design strategy has been chose, please validate this entry: ", template_information[template]["template_name"])
        
        # run one optimization step
        model_to_optimize.run(num_recycles=num_recycles, num_models=num_models, sample_models=sample_models,
             models=models, backprop=backprop, callback=callback)
        
        grad_t = copy.deepcopy(model_to_optimize.aux["grad"]['seq'])
        if template_information[template]["template_targeting_strategy"] == "detarget" and decay_scale is not None:
            grad_t = decay_scale * grad_t
        grads.append(grad_t)
        # float, not int: a fractional weight used to raise ValueError as a string
        # ("0.5") or silently truncate to 0 as a number, dropping the template
        try:
            base_ws.append(float(template_information[template]["loss_weighting"]))
        except (TypeError, ValueError):
            base_ws.append(1.0)

        # store gradient and loss
        template_information[template]["logging"]["loss_1"] = model_to_optimize.aux["log"]['loss']
        template_information[template]["logging"]["loss_2"] = model_to_optimize.aux['loss']
        
        template_information[template]["logging"]["pae_1"] = copy.deepcopy(model_to_optimize.aux["log"]['pae'])
        template_information[template]["logging"]["pae_2"] = copy.deepcopy(model_to_optimize.aux['losses']['pae'])

        template_information[template]["logging"]["i_pae_1"] = copy.deepcopy(model_to_optimize.aux["log"]['i_pae'])
        template_information[template]["logging"]["i_pae_2"] = copy.deepcopy(model_to_optimize.aux['losses']['i_pae'])

        template_information[template]["logging"]["con_1"] = copy.deepcopy(model_to_optimize.aux["log"]['con'])
        template_information[template]["logging"]["con_2"] = copy.deepcopy(model_to_optimize.aux['losses']['con'])

        template_information[template]["logging"]["i_con_1"] = copy.deepcopy(model_to_optimize.aux["log"]['i_con'])
        template_information[template]["logging"]["i_con_2"] = copy.deepcopy(model_to_optimize.aux['losses']['i_con'])

        template_information[template]["logging"]["plddt_1"] = copy.deepcopy(model_to_optimize.aux["log"]['plddt'])
        template_information[template]["logging"]["plddt_2"] = copy.deepcopy(model_to_optimize.aux['losses']['plddt'])

        template_information[template]["logging"]["ptm_1"] = copy.deepcopy(model_to_optimize.aux["log"]['ptm'])

        template_information[template]["logging"]["i_ptm_1"] = copy.deepcopy(model_to_optimize.aux["log"]['i_ptm'])
        # setting up the loss for the target plddt. 
        # template_information[template]["logging"]["plddt_target"] = copy.deepcopy(model_to_optimize.aux["plddt"][:-binder_length])

    # ----- aggregate the per-template gradients -----
    base_ws = np.asarray(base_ws, dtype=np.float32)
    losses_1 = np.asarray([template_information[t]["logging"]["loss_1"] for t in template_order],
                          dtype=np.float32)
    losses_2 = np.asarray([template_information[t]["logging"]["loss_2"] for t in template_order],
                          dtype=np.float32)

    # optional dynamic reweighting towards worst templates (uses logged loss_1)
    # dyn_beta=0 disables (weights remain base_ws). The softmax sums to 1, so enabling
    # this also scales the aggregate gradient down by roughly the template count.
    if dyn_beta is not None and float(dyn_beta) > 0.0 and len(losses_1) > 1:
        # shift improves numerical stability and makes beta scale more interpretable
        L_shift = losses_1 - np.mean(losses_1)
        dyn = softmax(float(dyn_beta) * L_shift)
        ws = base_ws * dyn
    else:
        ws = base_ws

    # avoid degenerate all-zero weights
    if not np.isfinite(ws).all() or float(np.sum(ws)) <= 0.0:
        ws = np.ones_like(base_ws, dtype=np.float32)

    running_gradient = None
    for grad_t, w in zip(grads, ws):
        if grad_rms_norm:
            # put every template on the same gradient scale, so loss_weighting rather
            # than raw gradient magnitude decides the balance between them
            grad_t = grad_t / (jnp.sqrt(jnp.mean(jnp.square(grad_t))) + grad_rms_eps)
        grad_t = grad_t * float(w)
        running_gradient = grad_t if running_gradient is None else running_gradient + grad_t

    current_sequence = copy.deepcopy(model_to_optimize.get_seq()[0])
    model_to_optimize.aux["grad"]['seq'] = running_gradient
    if model_to_optimize.opt["norm_seq_grad"]: model_to_optimize._norm_seq_grad()
    model_to_optimize._state, model_to_optimize.aux["grad"] = model_to_optimize._optimizer(model_to_optimize._state, model_to_optimize.aux["grad"], model_to_optimize._params)
    # apply gradients
    lr = model_to_optimize.opt["learning_rate"] * lr_scale
    model_to_optimize._params = jax.tree_util.tree_map(lambda x,g:x-lr*g, model_to_optimize._params, model_to_optimize.aux["grad"])
    
    # combine log/losses for tracking:
    # model_to_optimize.aux["log"]['loss'] += copy.deepcopy(np.sum([template_information[template]["logging"]["loss_1"] for template in template_information]))
    # model_to_optimize.aux['loss'] += copy.deepcopy(np.sum([template_information[template]["logging"]["loss_2"] for template in template_information]))
    # Tracking loss: weighted mean (default) or soft bottleneck (softmax over
    # per-template losses, approaching the worst template as alpha grows).
    # best_metric='loss' selects the saved design on this value.
    ws_sum = float(np.sum(ws)) if float(np.sum(ws)) > 0 else 1.0
    weighted_mean_loss1 = float(np.sum(ws * losses_1) / ws_sum)
    weighted_mean_loss2 = float(np.sum(ws * losses_2) / ws_sum)

    if bottleneck_alpha is not None and float(bottleneck_alpha) > 0.0 and len(losses_1) > 1:
        p = softmax(float(bottleneck_alpha) * (losses_1 - np.max(losses_1)))
        track_loss1 = float(np.sum(p * losses_1))
        track_loss2 = float(np.sum(p * losses_2))
    else:
        track_loss1 = weighted_mean_loss1
        track_loss2 = weighted_mean_loss2

    model_to_optimize.aux["log"]['loss'] = copy.deepcopy(track_loss1)
    model_to_optimize.aux["loss"] = copy.deepcopy(track_loss2)
    
    model_to_optimize.aux["log"]['pae'] = copy.deepcopy(np.mean([template_information[template]["logging"]["pae_1"] for template in template_information]))
    model_to_optimize.aux["losses"]['pae'] = copy.deepcopy(np.mean([template_information[template]["logging"]["pae_2"] for template in template_information]))

    model_to_optimize.aux["log"]['i_pae'] = copy.deepcopy(np.mean([template_information[template]["logging"]["i_pae_1"] for template in template_information]))
    model_to_optimize.aux["losses"]['i_pae'] = copy.deepcopy(np.mean([template_information[template]["logging"]["i_pae_2"] for template in template_information]))

    model_to_optimize.aux["log"]['con'] = copy.deepcopy(np.mean([template_information[template]["logging"]["con_1"] for template in template_information]))
    model_to_optimize.aux["losses"]['con'] = copy.deepcopy(np.mean([template_information[template]["logging"]["con_2"] for template in template_information]))

    model_to_optimize.aux["log"]['i_con'] = copy.deepcopy(np.mean([template_information[template]["logging"]["i_con_1"] for template in template_information]))
    model_to_optimize.aux["losses"]['i_con'] = copy.deepcopy(np.mean([template_information[template]["logging"]["i_con_2"] for template in template_information]))

    model_to_optimize.aux["log"]['plddt'] = copy.deepcopy(np.mean([template_information[template]["logging"]["plddt_1"] for template in template_information]))
    model_to_optimize.aux["losses"]['plddt'] = copy.deepcopy(np.mean([template_information[template]["logging"]["plddt_2"] for template in template_information]))

    model_to_optimize.aux["log"]['ptm'] = copy.deepcopy(np.mean([template_information[template]["logging"]["ptm_1"] for template in template_information]))

    model_to_optimize.aux["log"]['i_ptm'] = copy.deepcopy(np.mean([template_information[template]["logging"]["i_ptm_1"] for template in template_information]))

    # model_to_optimize.aux["log"]['sequences'] = current_sequence
    # tracking loss for the plddt of the target.
    # model_to_optimize.aux["log"]["plddt_target"] = copy.deepcopy(np.mean([template_information[template]["logging"]["plddt_target"] for template in template_information]))
    
    # save results
    model_to_optimize._save_results(save_best=save_best, verbose=verbose)
    # increment
    model_to_optimize._k += 1
   
   


def design_multi(model_to_optimize, template_information, advanced_settings, iters=100,
             soft=0.0, e_soft=None,
             temp=1.0, e_temp=None,
             hard=0.0, e_hard=None,
             step=1.0, e_step=None,
             dropout=True, opt=None, weights=None, 
             num_recycles=None, ramp_recycles=False, 
             num_models=None, sample_models=None, models=None,
             backprop=True, callback=None, save_best=False, verbose=1, binder_length=None, decay_iters=None, decay_scale=None):
    # update options/settings (if defined)

    # print(binder_length)

    model_to_optimize.set_opt(opt, dropout=dropout)
    model_to_optimize.set_weights(weights)

    m = {"soft":[soft,e_soft],"temp":[temp,e_temp],
         "hard":[hard,e_hard],"step":[step,e_step]}
    m = {k:[s,(s if e is None else e)] for k,(s,e) in m.items()}

    if ramp_recycles:
        if num_recycles is None:
            num_recycles = model_to_optimize.opt["num_recycles"]
        m["num_recycles"] = [0,num_recycles]
    
    for i in range(iters):
        for k,(s,e) in m.items():
            if k == "temp":
                model_to_optimize.set_opt({k:(e+(s-e)*(1-(i+1)/iters)**2)})
            else:
                v = (s+(e-s)*((i+1)/iters))
                if k == "step": step = v
                elif k == "num_recycles": num_recycles = round(v)
                else: model_to_optimize.set_opt({k:v})

        # decay learning rate based on temperature
        lr_scale = step * ((1 - model_to_optimize.opt["soft"]) + (model_to_optimize.opt["soft"] * model_to_optimize.opt["temp"]))   

        step_multi(model_to_optimize, template_information, advanced_settings=advanced_settings, lr_scale=lr_scale, num_recycles=num_recycles, num_models=num_models, sample_models=sample_models, models=models, backprop=backprop, callback=callback, save_best=save_best, verbose=verbose, binder_length=binder_length, iter_idx=i, decay_iters=decay_iters, decay_scale=decay_scale)


def design_mcmc_multi(model_to_optimize, template_information, advanced_settings, steps=500, half_life=200, T_init=0.01, mutation_rate=1, seq_logits=None, save_best=True, sample_models=True, verbose=1,  **kwargs):
   # gather settings
    model_flags = {k:kwargs.pop(k,None) for k in ["num_models","sample_models","models"]}

    # initialize with set numbers
    plddt, best_loss, current_loss = None, np.inf, np.inf 
    current_seq = (model_to_optimize._params["seq"] + model_to_optimize._inputs["bias"]).argmax(-1)
    if seq_logits is None: seq_logits = 0

    for template in template_information:  
        template_information[template]["MCMC_loss_value_logging"] = {0: np.inf}

    # run MCMC algorithm
    print("Running MCMC with simulated annealing...")
    for i in range(steps):

        # update temperature
        T = T_init * (np.exp(np.log(0.5) / half_life) ** i) 

        # mutate sequence
        if i == 0:
            mut_seq = current_seq
            last_accepted_index = 0
        else:
            mut_seq = model_to_optimize._mutate(
                seq=current_seq, 
                plddt=plddt,
                logits=seq_logits + model_to_optimize._inputs["bias"],
                mutation_rate=mutation_rate
            )
        # get loss for each of the templates that we have
        for template in template_information:
            
            # if template_information[template]["template_targeting_strategy"] != "target":
            #     continue
            
            # swap batch for current template:
            if template_information[template]["batch"] != None:
                model_to_optimize._inputs['batch'] = template_information[template]["batch"]
            # swap the hotspot too: prep_inputs() left opt["hotspot"] on whichever template was prepped last
            if template_information[template].get("hotspot") is None:
                model_to_optimize.opt.pop("hotspot", None)
            else:
                model_to_optimize.opt["hotspot"] = template_information[template]["hotspot"]

            model_nums = model_to_optimize._get_model_nums(**model_flags)
            aux = model_to_optimize.predict(seq=mut_seq, return_aux=True, verbose=False, model_nums=model_nums, **kwargs)
            loss = aux["log"]["loss"]

            template_information[template]["logging"]["loss_1"] = model_to_optimize.aux["log"]['loss']
            template_information[template]["logging"]["loss_2"] = model_to_optimize.aux['loss']
            
            template_information[template]["logging"]["pae_1"] = copy.deepcopy(model_to_optimize.aux["log"]['pae'])
            template_information[template]["logging"]["pae_2"] = copy.deepcopy(model_to_optimize.aux['losses']['pae'])

            template_information[template]["logging"]["i_pae_1"] = copy.deepcopy(model_to_optimize.aux["log"]['i_pae'])
            template_information[template]["logging"]["i_pae_2"] = copy.deepcopy(model_to_optimize.aux['losses']['i_pae'])

            template_information[template]["logging"]["con_1"] = copy.deepcopy(model_to_optimize.aux["log"]['con'])
            template_information[template]["logging"]["con_2"] = copy.deepcopy(model_to_optimize.aux['losses']['con'])

            template_information[template]["logging"]["i_con_1"] = copy.deepcopy(model_to_optimize.aux["log"]['i_con'])
            template_information[template]["logging"]["i_con_2"] = copy.deepcopy(model_to_optimize.aux['losses']['i_con'])

            template_information[template]["logging"]["plddt_1"] = copy.deepcopy(model_to_optimize.aux["log"]['plddt'])
            template_information[template]["logging"]["plddt_2"] = copy.deepcopy(model_to_optimize.aux['losses']['plddt'])

            template_information[template]["logging"]["ptm_1"] = copy.deepcopy(model_to_optimize.aux["log"]['ptm'])

            template_information[template]["logging"]["i_ptm_1"] = copy.deepcopy(model_to_optimize.aux["log"]['i_ptm'])


            if template_information[template]["template_targeting_strategy"] == "target":
                template_information[template]["MCMC_loss_value_logging"][i+1] = loss * int(template_information[template]["loss_weighting"])
            if template_information[template]["template_targeting_strategy"] == "detarget":
                # add a buffer for detargeting, as we want to allow for slight increase in loss 
                template_information[template]["MCMC_loss_value_logging"][i+1] = (loss - 0.1) * int(template_information[template]["loss_weighting"])

        total_loss_prev = sum(
            template_information[template]["MCMC_loss_value_logging"][last_accepted_index] 
            for template in template_information
        )
        total_loss_current = sum(
            template_information[template]["MCMC_loss_value_logging"][i+1]
            for template in template_information
        )

        loss_delta = total_loss_prev - total_loss_current

        # target_improved = (delta_target < 0)
        # if delta_detarget: detarget_tradeoff = (delta_detarget > 0 and delta_previous >= 0)
        # annealing_accept = (np.random.uniform() < np.exp(-delta_target / T))
        # # if annealing_accept: print(annealing_accept)

        # if accept_on_first or target_improved or detarget_tradeoff or annealing_accept:
        if loss_delta > 0 or i == 0:
            print("accepted MCMC:", loss)
            # accept
            last_accepted_index = i + 1
            (current_seq,current_loss) = (mut_seq,loss)
            
            plddt = aux["all"]["plddt"].mean(0)
            plddt = plddt[model_to_optimize._target_len:] if model_to_optimize.protocol == "binder" else plddt[:model_to_optimize._len]
            
            if loss < best_loss:
                (best_loss, model_to_optimize._k) = (loss, i)
                model_to_optimize.set_seq(seq=current_seq, bias=model_to_optimize._inputs["bias"])

                model_to_optimize.aux["log"]['loss'] = copy.deepcopy(np.mean([template_information[template]["logging"]["loss_1"] for template in template_information]))
                model_to_optimize.aux["loss"] = copy.deepcopy(np.mean([template_information[template]["logging"]["loss_2"] for template in template_information]))
                
                model_to_optimize.aux["log"]['pae'] = copy.deepcopy(np.mean([template_information[template]["logging"]["pae_1"] for template in template_information]))
                model_to_optimize.aux["losses"]['pae'] = copy.deepcopy(np.mean([template_information[template]["logging"]["pae_2"] for template in template_information]))

                model_to_optimize.aux["log"]['i_pae'] = copy.deepcopy(np.mean([template_information[template]["logging"]["i_pae_1"] for template in template_information]))
                model_to_optimize.aux["losses"]['i_pae'] = copy.deepcopy(np.mean([template_information[template]["logging"]["i_pae_2"] for template in template_information]))

                model_to_optimize.aux["log"]['con'] = copy.deepcopy(np.mean([template_information[template]["logging"]["con_1"] for template in template_information]))
                model_to_optimize.aux["losses"]['con'] = copy.deepcopy(np.mean([template_information[template]["logging"]["con_2"] for template in template_information]))

                model_to_optimize.aux["log"]['i_con'] = copy.deepcopy(np.mean([template_information[template]["logging"]["i_con_1"] for template in template_information]))
                model_to_optimize.aux["losses"]['i_con'] = copy.deepcopy(np.mean([template_information[template]["logging"]["i_con_2"] for template in template_information]))

                model_to_optimize.aux["log"]['plddt'] = copy.deepcopy(np.mean([template_information[template]["logging"]["plddt_1"] for template in template_information]))
                model_to_optimize.aux["losses"]['plddt'] = copy.deepcopy(np.mean([template_information[template]["logging"]["plddt_2"] for template in template_information]))

                model_to_optimize.aux["log"]['ptm'] = copy.deepcopy(np.mean([template_information[template]["logging"]["ptm_1"] for template in template_information]))

                model_to_optimize.aux["log"]['i_ptm'] = copy.deepcopy(np.mean([template_information[template]["logging"]["i_ptm_1"] for template in template_information]))


                model_to_optimize._save_results(save_best=save_best, verbose=True)

    return 0 


    
def design_semigreedy_multi(model_to_optimize, template_information, advanced_settings, iters=20, tries=8, dropout=False,
                        save_best=True, seq_logits=None, e_tries=None, binder_length=None, **kwargs):
    '''semigreedy search''' 
    if e_tries is None: e_tries = tries

    # get starting sequence
    if hasattr(model_to_optimize,"aux"):
        seq = model_to_optimize.aux["seq"]["logits"].argmax(-1)
    else:
        seq = (model_to_optimize._params["seq"] + model_to_optimize._inputs["bias"]).argmax(-1)


    # bias sampling towards the defined bias
    if seq_logits is None: seq_logits = 0
    
    model_flags = {k:kwargs.pop(k,None) for k in ["num_models","sample_models","models"]}
    verbose = kwargs.pop("verbose",1)

    # get current plddt
    aux = model_to_optimize.predict(seq, return_aux=True, verbose=False, **model_flags, **kwargs)
    plddt = model_to_optimize.aux["plddt"]
    plddt = plddt[model_to_optimize._target_len:] if model_to_optimize.protocol == "binder" else plddt[:model_to_optimize._len]

    # optimize!
    if verbose:
        print("Running semigreedy optimization...")
    
    for i in range(iters):
        buff = []
        model_nums = model_to_optimize._get_model_nums(**model_flags)
        num_tries = (tries+(e_tries-tries)*((i+1)/iters))
        for t in range(int(num_tries)):
            mut_seq = model_to_optimize._mutate(seq=seq, plddt=plddt, logits=seq_logits + model_to_optimize._inputs["bias"])

            for template in template_information:
                # empty & set template losses:
                template_information[template]["logging"] = {}

                # swap batch for current template:
                if template_information[template]["batch"] != None:
                    model_to_optimize._inputs['batch'] = template_information[template]["batch"]
                # swap the hotspot too: prep_inputs() left opt["hotspot"] on whichever template was prepped last
                if template_information[template].get("hotspot") is None:
                    model_to_optimize.opt.pop("hotspot", None)
                else:
                    model_to_optimize.opt["hotspot"] = template_information[template]["hotspot"]
                
                # swap weights for targeting/detargeting:
                if template_information[template]["template_targeting_strategy"] == "target":
                    model_to_optimize.opt["weights"].update(
                        {"pae":(advanced_settings["weights_pae_intra"]   ),
                        "plddt":(advanced_settings["weights_plddt"]      ),
                        "i_pae":(advanced_settings["weights_pae_inter"]  ),
                        "con":(advanced_settings["weights_con_intra"]    ),
                        "i_con":(advanced_settings["weights_con_inter"]  ),
                        "i_ptm":(advanced_settings["weights_iptm"]       ),
                        })
                    # print("targeting weights: ", model_to_optimize.opt["weights"])
                # scaling the weights according to detargeting scale
                elif template_information[template]["template_targeting_strategy"] == "detarget":
                    model_to_optimize.opt["weights"].update(
                        {"pae":(advanced_settings["weights_pae_intra"]   *(advanced_settings["detargeting_scaling_weights_pae_intra"])),
                        "plddt":(advanced_settings["weights_plddt"]      *(advanced_settings["detargeting_scaling_weights_plddt"])),
                        "i_pae":(advanced_settings["weights_pae_inter"]  *(advanced_settings["detargeting_scaling_weights_pae_inter"])),
                        "con":(advanced_settings["weights_con_intra"]    *(advanced_settings["detargeting_scaling_weights_con_intra"])),
                        "i_con":(advanced_settings["weights_con_inter"]  *(advanced_settings["detargeting_scaling_weights_con_inter"])),
                        "i_ptm":(advanced_settings["weights_iptm"]       *(advanced_settings["detargeting_scaling_weights_iptm"])),
            })
                    # print("detargeting weights: ", model_to_optimize.opt["weights"])
                else:
                    raise ValueError("No senseful design strategy has been chose, please validate this entry: ", template_information[template]["template_name"])
                

                # run one optimization step
                aux = model_to_optimize.predict(seq=mut_seq, return_aux=True, model_nums=model_nums, verbose=False, **kwargs)

                # store gradient and loss
                template_information[template]["logging"]["loss_1"] = model_to_optimize.aux["log"]['loss']
                template_information[template]["logging"]["loss_2"] = model_to_optimize.aux['loss']
                
                template_information[template]["logging"]["pae_1"] = copy.deepcopy(model_to_optimize.aux["log"]['pae'])
                template_information[template]["logging"]["pae_2"] = copy.deepcopy(model_to_optimize.aux['losses']['pae'])

                template_information[template]["logging"]["i_pae_1"] = copy.deepcopy(model_to_optimize.aux["log"]['i_pae'])
                template_information[template]["logging"]["i_pae_2"] = copy.deepcopy(model_to_optimize.aux['losses']['i_pae'])

                template_information[template]["logging"]["con_1"] = copy.deepcopy(model_to_optimize.aux["log"]['con'])
                template_information[template]["logging"]["con_2"] = copy.deepcopy(model_to_optimize.aux['losses']['con'])

                template_information[template]["logging"]["i_con_1"] = copy.deepcopy(model_to_optimize.aux["log"]['i_con'])
                template_information[template]["logging"]["i_con_2"] = copy.deepcopy(model_to_optimize.aux['losses']['i_con'])

                template_information[template]["logging"]["plddt_1"] = copy.deepcopy(model_to_optimize.aux["log"]['plddt'])
                template_information[template]["logging"]["plddt_2"] = copy.deepcopy(model_to_optimize.aux['losses']['plddt'])

                template_information[template]["logging"]["ptm_1"] = copy.deepcopy(model_to_optimize.aux["log"]['ptm'])

                template_information[template]["logging"]["i_ptm_1"] = copy.deepcopy(model_to_optimize.aux["log"]['i_ptm'])
                # setting up the loss for the target plddt. 
                # template_information[template]["logging"]["plddt_target"] = copy.deepcopy(model_to_optimize.aux["plddt"][:-binder_length])
            
            buff.append({"aux":aux, "seq":np.array(mut_seq)})

            # combine log/losses for tracking:
            # model_to_optimize.aux["log"]['loss'] += copy.deepcopy(np.sum([template_information[template]["logging"]["loss_1"] for template in template_information]))
            # model_to_optimize.aux['loss'] += copy.deepcopy(np.sum([template_information[template]["logging"]["loss_2"] for template in template_information]))

            model_to_optimize.aux["log"]['loss'] = copy.deepcopy(np.mean([template_information[template]["logging"]["loss_1"] for template in template_information]))
            model_to_optimize.aux['loss'] = copy.deepcopy(np.mean([template_information[template]["logging"]["loss_2"] for template in template_information]))

            model_to_optimize.aux["log"]['pae'] = copy.deepcopy(np.mean([template_information[template]["logging"]["pae_1"] for template in template_information]))
            model_to_optimize.aux["losses"]['pae'] = copy.deepcopy(np.mean([template_information[template]["logging"]["pae_2"] for template in template_information]))

            model_to_optimize.aux["log"]['i_pae'] = copy.deepcopy(np.mean([template_information[template]["logging"]["i_pae_1"] for template in template_information]))
            model_to_optimize.aux["losses"]['i_pae'] = copy.deepcopy(np.mean([template_information[template]["logging"]["i_pae_2"] for template in template_information]))

            model_to_optimize.aux["log"]['con'] = copy.deepcopy(np.mean([template_information[template]["logging"]["con_1"] for template in template_information]))
            model_to_optimize.aux["losses"]['con'] = copy.deepcopy(np.mean([template_information[template]["logging"]["con_2"] for template in template_information]))

            model_to_optimize.aux["log"]['i_con'] = copy.deepcopy(np.mean([template_information[template]["logging"]["i_con_1"] for template in template_information]))
            model_to_optimize.aux["losses"]['i_con'] = copy.deepcopy(np.mean([template_information[template]["logging"]["i_con_2"] for template in template_information]))

            model_to_optimize.aux["log"]['plddt'] = copy.deepcopy(np.mean([template_information[template]["logging"]["plddt_1"] for template in template_information]))
            model_to_optimize.aux["losses"]['plddt'] = copy.deepcopy(np.mean([template_information[template]["logging"]["plddt_2"] for template in template_information]))

            model_to_optimize.aux["log"]['ptm'] = copy.deepcopy(np.mean([template_information[template]["logging"]["ptm_1"] for template in template_information]))

            model_to_optimize.aux["log"]['i_ptm'] = copy.deepcopy(np.mean([template_information[template]["logging"]["i_ptm_1"] for template in template_information]))
            
            # tracking loss for the plddt of the target. 
            # model_to_optimize.aux["log"]["plddt_target"] = copy.deepcopy(np.mean([template_information[template]["logging"]["plddt_target"] for template in template_information]))


        # accept best
        losses = [x["aux"]["loss"] for x in buff]
        best = buff[np.argmin(losses)]
        model_to_optimize.aux, seq = best["aux"], jnp.array(best["seq"])
        model_to_optimize.set_seq(seq=seq, bias=model_to_optimize._inputs["bias"])
        model_to_optimize._save_results(save_best=save_best, verbose=verbose)
        # update plddt
        plddt = best["aux"]["plddt"]
        plddt = plddt[model_to_optimize._target_len:] if model_to_optimize.protocol == "binder" else plddt[:model_to_optimize._len]
        model_to_optimize._k += 1
    

def process_trajectory_complex(
    binder_sequence,
    design_name,         # Replaces mpnn_design_name
    target_pdb,          # Path to the target PDB file
    chain,               # Chain identifier
    length,              # Binder length
    trajectory_dir,      # Directory where trajectory PDBs will be saved
    prediction_models,
    failure_csv,         # Non-default argument comes before default arguments
    appendix="",   
    advanced_settings={} # Default argument placed at the end
):
    """
    Predicts the trajectory complex using a binder sequence and relaxes the generated models.
    
    Parameters:
      prediction_model: Base prediction model (used for metrics reference)
      binder_sequence: Amino acid sequence for the binder (will be cleaned)
      design_name: Name/identifier for the design (used in file naming)
      target_pdb: File path to the target PDB structure
      chain: Chain identifier for the target
      length: Length of the binder sequence
      trajectory_dir: Directory path to save trajectory PDB files
      relaxed_dir: Directory path to save relaxed PDB files
      prediction_models: Iterable of model numbers to predict/relax
      failure_csv: Path to the failure CSV
      appendix: Additional string to append to filenames
      advanced_settings: Dictionary containing advanced settings, e.g.:
           - "num_recycles_validation": int
           - "af_params_dir": str (path to AF parameters)
           - "rm_template_seq_predict": bool
           - "rm_template_sc_predict": bool
           - "multimer_validation": bool
      seed: Optional random seed for reproducibility
      
    Returns:
      prediction_stats: Dictionary containing prediction metrics for each model.
    """
    # Initialize the trajectory prediction model using advanced settings.
    trajectory_prediction_model = mk_afdesign_model(
        protocol="binder",
        num_recycles=advanced_settings.get("num_recycles_design", 1),
        data_dir=advanced_settings.get("af_params_dir", ""),
        use_multimer=advanced_settings.get("use_multimer_design", True)
    )
    
    # Prepare the inputs using the target PDB and chain.
    trajectory_prediction_model.prep_inputs(
        pdb_filename=target_pdb,
        chain=chain,
        binder_len=length,
        rm_target_seq=advanced_settings.get("rm_template_seq_design", False),
        rm_target_sc=advanced_settings.get("rm_template_sc_design", False)
    )
    
    prediction_stats = {}
    
    # Clean binder sequence: keep only uppercase letters.
    binder_sequence = re.sub("[^A-Z]", "", binder_sequence.upper())
    
    # Loop over each model number in the prediction_models list.
    for model_num in prediction_models:
        # Construct the file path for the trajectory PDB.
        trajectory_pdb = os.path.join(trajectory_dir, f"{design_name}{appendix}.pdb")
        
        # Predict the model only if the file does not already exist.
        if not os.path.exists(trajectory_pdb):
            trajectory_prediction_model.predict(
                seq=binder_sequence,
                models=[model_num],
                num_recycles=advanced_settings.get("num_recycles_design", 1),
                verbose=False
            )
            # Save the predicted structure.
            trajectory_prediction_model.save_pdb(trajectory_pdb)
            
            # Copy the prediction metrics from the model's auxiliary log.
            prediction_metrics = trajectory_prediction_model.aux["log"].copy()  # Expected keys: plddt, ptm, i_ptm, pae, i_pae
            
            # Store rounded statistics.
            stats = {
                'pLDDT': round(prediction_metrics.get('plddt', 0), 2), 
                'pTM': round(prediction_metrics.get('ptm', 0), 2), 
                'i_pTM': round(prediction_metrics.get('i_ptm', 0), 2), 
                'pAE': round(prediction_metrics.get('pae', 0), 2), 
                'i_pAE': round(prediction_metrics.get('i_pae', 0), 2),
                # 'plddt_target': round(prediction_metrics.get('plddt_target', 0), 2)
            }
            prediction_stats[model_num+1] = stats
            print(f"Model {model_num+1} stats:", stats)
        else:
            print(f"Trajectory file {trajectory_pdb} already exists. Skipping prediction.")
        
        ca_clashes = calculate_clash_score(trajectory_pdb, 2.5, only_ca=True)
        if ca_clashes > 0:
            trajectory_prediction_model.aux['log']['terminate'] = "Clashing"
            update_failures(failure_csv, 'Trajectory_Clashes')
            print("Severe clashes detected, skipping analysis and MPNN optimisation")
            print("")
        else:
            # check if low quality prediction
            if stats['pLDDT'] < 0.7:
                trajectory_prediction_model.aux['log']['terminate'] = "LowConfidence"
                update_failures(failure_csv, 'Trajectory_final_pLDDT')
                print("Trajectory starting confidence low, skipping analysis and MPNN optimisation")
                print("")
            else:
                # does it have enough contacts to consider?
                binder_contacts = hotspot_residues(trajectory_pdb)
                binder_contacts_n = len(binder_contacts.items())

                # if less than 7 contacts then protein is floating above and is not binder
                if binder_contacts_n < 7:
                    trajectory_prediction_model.aux['log']['terminate'] = "LowConfidence"
                    update_failures(failure_csv, 'Trajectory_Contacts')
                    print("Too few contacts at the interface, skipping analysis and MPNN optimisation")
                    print("")
                else:
                    # phew, trajectory is okay! We can continue
                    trajectory_prediction_model.aux["log"]["terminate"] = ""
                    print(f"Trajectory against {appendix} successful, final pLDDT: "+str(stats['pLDDT']))    
    
    return trajectory_prediction_model, trajectory_pdb

# Define radius of gyration loss for colabdesign
def add_rg_loss(self, weight=0.1):
    '''add radius of gyration loss'''
    def loss_fn(inputs, outputs):
        xyz = outputs["structure_module"]
        ca = xyz["final_atom_positions"][:,residue_constants.atom_order["CA"]]
        ca = ca[-self._binder_len:]
        rg = jnp.sqrt(jnp.square(ca - ca.mean(0)).sum(-1).mean() + 1e-8)
        rg_th = 2.38 * ca.shape[0] ** 0.365

        rg = jax.nn.elu(rg - rg_th)
        return {"rg":rg}

    self._callbacks["model"]["loss"].append(loss_fn)
    self.opt["weights"]["rg"] = weight
# Define interface pTM loss for colabdesign
def add_i_ptm_loss(self, weight=0.1):
    def loss_iptm(inputs, outputs):
        p = 1 - get_ptm(inputs, outputs, interface=True)
        i_ptm = mask_loss(p)
        return {"i_ptm": i_ptm}
    
    self._callbacks["model"]["loss"].append(loss_iptm)
    self.opt["weights"]["i_ptm"] = weight

# add target plddt loss
def add_target_plddt_loss(self, weight: float = 0.0) -> None:
    """Optimise predicted confidence over target residues."""
    def _per_residue_plddt(outputs) -> jnp.ndarray:
        logits = outputs["predicted_lddt"]["logits"]
        num_bins = logits.shape[-1]
        centres = jnp.arange(0.5 / num_bins, 1.0, 1.0 / num_bins)
        probs = jax.nn.softmax(logits, axis=-1)
        return jnp.sum(probs * centres, axis=-1)

    def loss_fn(inputs, outputs):
        plddt = _per_residue_plddt(outputs)
        target_plddt = plddt[:self._target_len]

        mean_val = target_plddt.mean()
        loss_val = 1.0 - mean_val
        return {"plddt_target": loss_val}

    self._callbacks["model"]["loss"].append(loss_fn)
    self.opt["weights"]["plddt_target"] = weight

# add helicity loss
def add_helix_loss(self, weight=0):
    def binder_helicity(inputs, outputs):  
      if "offset" in inputs:
        offset = inputs["offset"]
      else:
        idx = inputs["residue_index"].flatten()
        offset = idx[:,None] - idx[None,:]

      # define distogram
      dgram = outputs["distogram"]["logits"]
      dgram_bins = get_dgram_bins(outputs)
      mask_2d = np.outer(np.append(np.zeros(self._target_len), np.ones(self._binder_len)), np.append(np.zeros(self._target_len), np.ones(self._binder_len)))

      x = _get_con_loss(dgram, dgram_bins, cutoff=6.0, binary=True)
      if offset is None:
        if mask_2d is None:
          helix_loss = jnp.diagonal(x,3).mean()
        else:
          helix_loss = jnp.diagonal(x * mask_2d,3).sum() + (jnp.diagonal(mask_2d,3).sum() + 1e-8)
      else:
        mask = offset == 3
        if mask_2d is not None:
          mask = jnp.where(mask_2d,mask,0)
        helix_loss = jnp.where(mask,x,0.0).sum() / (mask.sum() + 1e-8)

      return {"helix":helix_loss}
    self._callbacks["model"]["loss"].append(binder_helicity)
    self.opt["weights"]["helix"] = weight

# Get pLDDT of best model
def get_best_plddt(af_model, length):
    return round(np.mean(af_model._tmp["best"]["aux"]["plddt"][-length:]),2)


# run MPNN to generate sequences for binders
def mpnn_gen_sequence(trajectory_pdb, binder_chain, trajectory_interface_residues, advanced_settings):
    # clear GPU memory
    clear_mem()

    # initialise MPNN model
    mpnn_model = mk_mpnn_model(backbone_noise=advanced_settings["backbone_noise"], model_name=advanced_settings["model_path"], weights=advanced_settings["mpnn_weights"])

    # check whether keep the interface generated by the trajectory or whether to redesign with MPNN
    design_chains = 'A,' + binder_chain

    if advanced_settings["mpnn_fix_interface"]:
        fixed_positions = 'A,' + trajectory_interface_residues
        print("Fixing interface residues on Binder: "+trajectory_interface_residues)
    else:
        fixed_positions = 'A'

    # prepare inputs for MPNN
    mpnn_model.prep_inputs(pdb_filename=trajectory_pdb, chain=design_chains, fix_pos=fixed_positions, rm_aa=advanced_settings["omit_AAs"])

    # sample MPNN sequences in parallel
    mpnn_sequences = mpnn_model.sample_parallel(temperature=advanced_settings["sampling_temp"], num=advanced_settings["num_seqs"], batch=advanced_settings["sample_seq_parallel"])

    return mpnn_sequences

# run prediction for binder with masked template target
def masked_binder_predict_multi(trajectory_dict, binder_sequence, mpnn_design_name, prediction_models, advanced_settings, filters, design_paths):
    prediction_stats = {}

    # clean sequence
    binder_sequence = re.sub("[^A-Z]", "", binder_sequence.upper())

    # reset filtering conditionals
    pass_af2_filters = True
    filter_failures = {}

    # start prediction per AF2 model, 2 are used by default due to masked templates
    for model_num in prediction_models:
        # check to make sure prediction does not exist already
        complex_pdb = os.path.join(design_paths["MPNN"], f"{mpnn_design_name}_model{model_num+1}{trajectory_dict['template_name']}.pdb")
        if not os.path.exists(complex_pdb):
            # predict model
            trajectory_dict["complex_prediction_model"].predict(seq=binder_sequence, models=[model_num], num_recycles=advanced_settings["num_recycles_validation"], verbose=False)
            trajectory_dict["complex_prediction_model"].save_pdb(complex_pdb)
            prediction_metrics = copy.deepcopy(trajectory_dict["complex_prediction_model"].aux["log"]) # contains plddt, ptm, i_ptm, pae, i_pae

            # extract the statistics for the model
            stats = {
                'pLDDT': round(prediction_metrics['plddt'], 2), 
                'pTM': round(prediction_metrics['ptm'], 2), 
                'i_pTM': round(prediction_metrics['i_ptm'], 2), 
                'pAE': round(prediction_metrics['pae'], 2), 
                'i_pAE': round(prediction_metrics['i_pae'], 2)
            }
            prediction_stats[model_num+1] = stats

            ## get filter conditions from the file

            # List of filter conditions and corresponding keys
            filter_conditions = [
                (f"{model_num+1}_pLDDT", 'plddt', '>='),
                (f"{model_num+1}_pTM", 'ptm', '>='),
                (f"{model_num+1}_i_pTM", 'i_ptm', '>='),
                (f"{model_num+1}_pAE", 'pae', '<='),
                (f"{model_num+1}_i_pAE", 'i_pae', '<='),
            ]

            # perform initial AF2 values filtering to determine whether to skip relaxation and interface scoring
            for filter_name, metric_key, _ in filter_conditions:
                threshold = filters.get(filter_name, {}).get("threshold")
                comparison = ">=" if filters.get(filter_name, {}).get("higher") else "<="
                if threshold is not None:
                    if comparison == '>=' and prediction_metrics[metric_key] < threshold:
                        pass_af2_filters = False
                        filter_failures[filter_name] = filter_failures.get(filter_name, 0) + 1
                    elif comparison == '<=' and prediction_metrics[metric_key] > threshold:
                        pass_af2_filters = False
                        filter_failures[filter_name] = filter_failures.get(filter_name, 0) + 1

            if not pass_af2_filters:
                break

    # Update the CSV file with the failure counts
    if filter_failures:
        update_failures(trajectory_dict["failure_csv"], filter_failures)

    # AF2 filters passed, contuing with relaxation
    for model_num in prediction_models:
        complex_pdb = os.path.join(design_paths["MPNN"], f"{mpnn_design_name}_model{model_num+1}{trajectory_dict['template_name']}.pdb")
        if pass_af2_filters:
            mpnn_relaxed = os.path.join(design_paths["MPNN/Relaxed"], f"{mpnn_design_name}_model{model_num+1}{trajectory_dict['template_name']}.pdb")
            pr_relax(complex_pdb, mpnn_relaxed)
        else:
            if os.path.exists(complex_pdb):
                os.remove(complex_pdb)

    return prediction_stats, pass_af2_filters


# run prediction for binder with masked template target
def masked_binder_predict_multi_adjusted(trajectory_dict, binder_sequence, mpnn_design_name, prediction_models, advanced_settings, filters, design_paths, target_statistics):
    prediction_stats = {}
    targeting_strategy = trajectory_dict["template_targeting_strategy"]

    # clean sequence
    binder_sequence = re.sub("[^A-Z]", "", binder_sequence.upper())

    # reset filtering conditionals
    pass_af2_filters = True
    filter_failures = {}

    # start prediction per AF2 model, 2 are used by default due to masked templates
    for model_num in prediction_models:
        # check to make sure prediction does not exist already
        complex_pdb = os.path.join(design_paths["MPNN"], f"{mpnn_design_name}_model{model_num+1}{trajectory_dict['template_name']}.pdb")
        if not os.path.exists(complex_pdb):
            # predict model
            trajectory_dict["complex_prediction_model"].predict(seq=binder_sequence, models=[model_num], num_recycles=advanced_settings["num_recycles_validation"], verbose=False)
            trajectory_dict["complex_prediction_model"].save_pdb(complex_pdb)
            prediction_metrics = copy.deepcopy(trajectory_dict["complex_prediction_model"].aux["log"]) # contains plddt, ptm, i_ptm, pae, i_pae

            # extract the statistics for the model
            stats = {
                'pLDDT': round(prediction_metrics['plddt'], 2), 
                'pTM': round(prediction_metrics['ptm'], 2), 
                'i_pTM': round(prediction_metrics['i_ptm'], 2), 
                'pAE': round(prediction_metrics['pae'], 2), 
                'i_pAE': round(prediction_metrics['i_pae'], 2)
            }
            prediction_stats[model_num+1] = stats

            ## get filter conditions from the file

            # List of filter conditions and corresponding keys
            filter_conditions = [
                (f"{model_num+1}_pLDDT", 'plddt', '>='),
                (f"{model_num+1}_pTM", 'ptm', '>='),
                (f"{model_num+1}_i_pTM", 'i_ptm', '>='),
                (f"{model_num+1}_pAE", 'pae', '<='),
                (f"{model_num+1}_i_pAE", 'i_pae', '<='),
            ]

            if targeting_strategy != "detarget":
                # perform initial AF2 values filtering to determine whether to skip relaxation and interface scoring
                for filter_name, metric_key, _ in filter_conditions:
                    threshold = filters.get(filter_name, {}).get("threshold")
                    comparison = ">=" if filters.get(filter_name, {}).get("higher") else "<="
                    if threshold is not None:
                        if comparison == '>=' and prediction_metrics[metric_key] < threshold:
                            pass_af2_filters = False
                            filter_failures[filter_name] = filter_failures.get(filter_name, 0) + 1
                        elif comparison == '<=' and prediction_metrics[metric_key] > threshold:
                            pass_af2_filters = False
                            filter_failures[filter_name] = filter_failures.get(filter_name, 0) + 1

            else: 
                for filter_name, metric_key, _ in filter_conditions:
                    if "i_pTM" in filter_name:
                        threshold = filters.get(filter_name, {}).get("threshold")
                        threshold = target_statistics["i_pTM"] - threshold
                        comparison = ">=" if filters.get(filter_name, {}).get("higher") else "<="
                        if threshold is not None:
                            if comparison == '>=' and prediction_metrics[metric_key] < threshold:
                                pass_af2_filters = False
                                filter_failures[filter_name] = filter_failures.get(filter_name, 0) + 1
                            elif comparison == '<=' and prediction_metrics[metric_key] > threshold:
                                pass_af2_filters = False
                                filter_failures[filter_name] = filter_failures.get(filter_name, 0) + 1
                    
                    else: 
                        threshold = filters.get(filter_name, {}).get("threshold")
                        comparison = ">=" if filters.get(filter_name, {}).get("higher") else "<="
                        if threshold is not None:
                            if comparison == '>=' and prediction_metrics[metric_key] < threshold:
                                pass_af2_filters = False
                                filter_failures[filter_name] = filter_failures.get(filter_name, 0) + 1
                            elif comparison == '<=' and prediction_metrics[metric_key] > threshold:
                                pass_af2_filters = False
                                filter_failures[filter_name] = filter_failures.get(filter_name, 0) + 1


            if not pass_af2_filters:
                break

    # Update the CSV file with the failure counts
    if filter_failures:
        update_failures(trajectory_dict["failure_csv"], filter_failures)

    # AF2 filters passed, contuing with relaxation
    for model_num in prediction_models:
        complex_pdb = os.path.join(design_paths["MPNN"], f"{mpnn_design_name}_model{model_num+1}{trajectory_dict['template_name']}.pdb")
        if pass_af2_filters:
            mpnn_relaxed = os.path.join(design_paths["MPNN/Relaxed"], f"{mpnn_design_name}_model{model_num+1}{trajectory_dict['template_name']}.pdb")
            pr_relax(complex_pdb, mpnn_relaxed)
        else:
            if os.path.exists(complex_pdb):
                os.remove(complex_pdb)

    return prediction_stats, pass_af2_filters

# run prediction for binder alone
def predict_binder_alone(prediction_model, binder_sequence, mpnn_design_name, length, trajectory_pdb, binder_chain, prediction_models, advanced_settings, design_paths, seed=None):
    binder_stats = {}

    # prepare sequence for prediction
    binder_sequence = re.sub("[^A-Z]", "", binder_sequence.upper())
    prediction_model.set_seq(binder_sequence)

    # predict each model separately
    for model_num in prediction_models:
        # check to make sure prediction does not exist already
        binder_alone_pdb = os.path.join(design_paths["MPNN/Binder"], f"{mpnn_design_name}_model{model_num+1}.pdb")
        if not os.path.exists(binder_alone_pdb):
            # predict model
            prediction_model.predict(models=[model_num], num_recycles=advanced_settings["num_recycles_validation"], verbose=False)
            prediction_model.save_pdb(binder_alone_pdb)
            prediction_metrics = copy_dict(prediction_model.aux["log"]) # contains plddt, ptm, pae

            # align binder model to trajectory binder
            align_pdbs(trajectory_pdb, binder_alone_pdb, binder_chain, "A")

            # extract the statistics for the model
            stats = {
                'pLDDT': round(prediction_metrics['plddt'], 2), 
                'pTM': round(prediction_metrics['ptm'], 2), 
                'pAE': round(prediction_metrics['pae'], 2)
            }
            binder_stats[model_num+1] = stats

    return binder_stats



# plot design trajectory losses
def plot_trajectory(af_model, design_name, design_paths):
    metrics_to_plot = ['loss', 'plddt', 'ptm', 'i_ptm', 'con', 'i_con', 'pae', 'i_pae', 'rg', 'mpnn', 'i_ptm_spread']
    colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k']

    metrics_dict = {}

    for index, metric in enumerate(metrics_to_plot):
        if metric in af_model.aux["log"]:
            # Create a new figure for each metric
            plt.figure()

            loss = af_model.get_loss(metric)
            metrics_dict[metric] = loss.tolist() if hasattr(loss, 'tolist') else list(loss)
            # Create an x axis for iterations
            iterations = range(1, len(loss) + 1)

            plt.plot(iterations, loss, label=f'{metric}', color=colors[index % len(colors)])

            # Add labels and a legend
            plt.xlabel('Iterations')
            plt.ylabel(metric)
            plt.title(design_name)
            plt.legend()
            plt.grid(True)

            # Save the plot
            plt.savefig(os.path.join(design_paths["Trajectory/Plots"], design_name+"_"+metric+".png"), dpi=150)
            
            # Close the figure
            plt.close()

    metrics_file_path = os.path.join(design_paths["Trajectory/Plots"], f"{design_name}_metrics.json")
    with open(metrics_file_path, 'w') as f:
        json.dump(metrics_dict, f, indent=4)
