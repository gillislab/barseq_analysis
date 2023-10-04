

source("dataset.R")


main = function() {
    sce = load_glu()
    sce$x = fct_rev(make_bins(sce$CCF_streamlines_x, 50))
    sce$y = fct_rev(make_bins(sce$CCF_streamlines_y, 50))
    sce$h2 = factor(sce$h2, ctx_subclass())
    sce$sample = sce$study_id
    
    # plot flatmaps / depth
    plot_ct_position(sce)
    
    # test for differences in composition between control and enucleated samples
    plot_differential_composition(sce)
}

load_glu = function() {
    labels = load_labels()
    result = lapply(set_names(sample_names()), function(my_sample) {
        sce = load_sample(my_sample)
        my_labels = labels[match(colnames(sce), labels$cell_id),]
        colData(sce) = cbind(colData(sce), select(my_labels, -cell_id))
        sce = sce[, sce$h1 == "Glutamatergic" & sce$h2 %in% ctx_subclass()]
        return(sce)
    })
    result = MetaNeighbor::mergeSCE(result)
    return(result)
}

plot_ct_position = function(sce) {
    plot_depth(sce)
    plot_flatmap(sce)
}

make_bins = function(x, n_bins) {
    x = cut(x, breaks = seq(min(x,na.rm = TRUE),max(x, na.rm = TRUE),length.out = n_bins))
}

plot_flatmap = function(sce) {
    h3_to_h2 = deframe(distinct(tibble(sce$h3_name, sce$h2)))

    h2_flatmap = compute_flatmap(sce, "h2")
    h3_flatmap = compute_flatmap(sce, "h3_name")
    h3_flatmap$supergroup = h3_to_h2[h3_flatmap$group]
    pdf("fig/flatmap.pdf")
    print(plot_flatmap_by_group(h2_flatmap, n_threshold=5))
    p = plot_flatmap_by_group(h3_flatmap, n_threshold = 2) +
        facet_wrap(~ supergroup + group) +
        theme(strip.text.x = element_text(size = 6))
    print(p)
    for (my_h2 in levels(sce$h2)) {
        print(plot_flatmap_by_group(filter(h3_flatmap,supergroup==my_h2), n_threshold=5))
    }
    dev.off()
    
    samples = sort(unique(sce$sample))
    h2_flatmap = bind_rows(lapply(set_names(samples), function(my_sample) {
        compute_flatmap(filter(as_tibble(colData(sce)), sample == my_sample), "h2")
    }), .id="sample")
    h3_flatmap = bind_rows(lapply(set_names(samples), function(my_sample) {
        fm = compute_flatmap(filter(as_tibble(colData(sce)), sample == my_sample), "h3_name")
        fm$supergroup = h3_to_h2[fm$group]
        fm
    }), .id="sample")
    pdf("fig/flatmap_by_sample.pdf")
    print(plot_flatmap_by_group(h2_flatmap, n_threshold=5) + facet_grid(group ~ sample))
    for (my_h2 in levels(sce$h2)) {
        print(plot_flatmap_by_group(filter(h3_flatmap,supergroup==my_h2), n_threshold=5) + facet_grid(group ~ sample))
    }
    dev.off()
}

compute_flatmap = function(sce, group_by, supergroup_by=NULL) {
    result = tibble(x=sce$x, y=sce$y, group = sce[[group_by]],
                    supergroup = if(is.null(supergroup_by)) {"all"} else {sce[[supergroup_by]]}) %>%
        filter(!is.na(x) & ! is.na(y)) %>%
        group_by(x,y,group,supergroup) %>%
        tally() %>%
        ungroup() %>%
        complete(nesting(group,supergroup), nesting(x,y), fill=list("n"=0)) %>%
        group_by(x,y,supergroup) %>%
        mutate(n_total = sum(n), f=n/n_total) %>%
        group_by(group, supergroup) %>%
        mutate(f01 = f / max(f)) %>%
        ungroup()
    return(result)
}

plot_flatmap_by_group= function(to_plot, n_threshold=0) {
    p = to_plot %>%
        mutate(f01 = ifelse(n_total >= n_threshold, f01, 0)) %>%
        ggplot(aes(x,y,fill=f01)) +
        geom_raster(show.legend = FALSE) +
        facet_wrap(~ group) +
        theme_void() +
        scale_fill_gradient(low = gray(0.95), high = "darkblue", limits=c(0,NA)) +
        coord_equal()
    return(p)
}

plot_depth = function(sce) {
    to_plot = tibble(depth=sce$CCF_streamlines_z, sample = sce$sample, h2 = sce$h2, h3 = sce$h3_name)
    p = to_plot %>%
        ggplot(aes(x=h2,y=-depth,fill = h2)) +
        geom_violin(show.legend = FALSE, scale = "width") +
        theme_classic(base_size = 16) +
        labs(x=NULL,y="Depth") +
        theme(axis.text.x = element_text(angle=45, hjust=1))
    ggsave("fig/depth_h2.pdf",p)
    p = to_plot %>%
        ggplot(aes(x=h2,y=-depth,fill = h2)) +
        geom_violin(show.legend = FALSE, scale = "width") +
        facet_wrap(~sample) +
        theme_classic(base_size = 12) +
        labs(x=NULL,y="Depth") +
        theme(axis.text.x = element_text(angle=45, hjust=1))
    ggsave("fig/depth_h2_by_sample.pdf",p)
    p = to_plot %>%
        ggplot(aes(x=sample,y=-depth,fill = h2)) +
        geom_violin(show.legend = FALSE, scale = "width") +
        facet_wrap(~h2) +
        theme_classic(base_size = 12) +
        labs(x=NULL,y="Depth") +
        theme(axis.text.x = element_text(angle=45, hjust=1))
    ggsave("fig/depth_sample_by_h2.pdf",p)
    p = to_plot %>%
        ggplot(aes(x=h3,y=-depth,fill = h2)) +
        geom_violin(show.legend = FALSE, scale = "width") +
        facet_wrap(~h2, scales="free_x") +
        theme_classic(base_size = 12) +
        labs(x=NULL,y="Depth") +
        theme(axis.text.x = element_text(angle=45, hjust=1))
    ggsave("fig/depth_h3.pdf",p)
}

plot_differential_compositon = function(sce) {
    my_data = tibble(sample = sce$sample, h2 = sce$h2, h3 = sce$h3_name, ccf1 = sce$CCFparentname, ccf2 = sce$CCFname) %>%
        filter(ccf1 %in% ctx_roi())

    ## H2 level
    h2 = aggregate_ct_ccf(my_data, "h2", "ccf1", min_min = 10)
    h2_aov = lapply(set_names(sort(unique(h2$h2))), function(my_ct) {lm(log2(prop) ~ litter*ccf1 + ccf1:condition, data = filter(h2, h2 == my_ct))})
    h2_varexp = summarize_varexp(h2_aov)
    h2_summary = summarize_results(h2_aov)
    h2_mean_prop = h2 %>%
        group_by(h2, ccf1) %>%
        summarize(mean_prop = mean(prop))
    qplot(h2_summary$pval, boundary=0) + theme_classic(base_size = 16) + labs(x="P-value (H2-level ANOVA)")
    ggsave("fig/diff_h2_pval.pdf")

    # MA plot
    h2_summary %>%
        filter(grepl("condition",variable)) %>%
        mutate(ccf1 = str_extract(variable, "^ccf1([^:]+)", group = 1)) %>%
        mutate(var = paste(cell_type, ccf1, sep = "|"), h2=cell_type) %>%
        left_join(h2_mean_prop) %>%
        mutate(fdr_bin = cut(sign(Estimate)*fdr, breaks = c(-1,-0.05,0,0.05,1))) %>%
        ggplot(aes(x=mean_prop,y=Estimate, label=var, color = fdr_bin)) +
        geom_point(size=0.5, show.legend = FALSE) +
        geom_hline(yintercept = 0, linetype="dashed") + 
        ggrepel::geom_text_repel(data= . %>% filter(fdr<0.05), size=5, show.legend = FALSE) +
        theme_classic(base_size = 16) +
        labs(x = "Average proportion", y="log2(FC)", col="FDR") +
        scale_x_log10() +
        scale_color_manual(values = c("gray70","blue3","brown1","gray70"))
    ggsave("fig/diff_2_ma.pdf")
    
    h2_significant = h2_summary %>%
        filter(grepl("condition",variable) & fdr<0.1) %>%
        mutate(ccf = str_extract(variable, "^ccf1([^:]+)", group = 1)) %>%
        mutate(var = paste(cell_type, ccf, sep = "|")) %>%
        mutate(var = fct_reorder(var, fdr))
    h2 %>%
        mutate(var = paste(h2, ccf1, sep="|")) %>%
        filter(var %in% filter(h2_significant, Estimate>0)$var) %>%
        mutate(var = factor(var, levels(h2_significant$var))) %>%
        ggplot(aes(x=litter,y=prop,fill=condition)) +
        geom_col(position = "dodge") +
        facet_wrap(~ var) +
        theme_classic(base_size=12) +
        labs(x = "Litter", y="Proportion of cortical excitatory neurons",title="Enriched in ENU")
    h2 %>%
        mutate(var = paste(h2, ccf1, sep="|")) %>%
        filter(var %in% filter(h2_significant, Estimate<0)$var) %>%
        mutate(var = factor(var, levels(h2_significant$var))) %>%
        ggplot(aes(x=litter,y=prop,fill=condition)) +
        geom_col(position = "dodge") +
        facet_wrap(~ var) +
        theme_classic(base_size=12) +
        labs(x = "Litter", y="Proportion of cortical excitatory neurons",title="Enriched in CTL")
    
    ## H3 level - fold change as effect size
    h3_to_h2 = deframe(distinct(select(my_data, h3, h2)))
    h3 = aggregate_ct_ccf(my_data, "h3", "ccf1", min_min = 10) %>%
        group_by(h3) %>%
        filter(length(unique(ccf1))>1)
    # try qlogis, log, log2 -> roughly equivalent, log2 more interpretable
    h3_aov = lapply(set_names(sort(unique(h3$h3))), function(my_ct) {lm(log2(prop+1e-5) ~ litter*ccf1 + ccf1:condition, data = filter(h3, h3 == my_ct))})
    h3_varexp = summarize_varexp(h3_aov)
    h3_summary = summarize_results(h3_aov)
    qplot(h3_summary$pval, boundary=0) + theme_classic(base_size = 16) + labs(x="P-value (H3-level ANOVA)")
    ggsave("fig/diff_h3_pval.pdf")
    
    # basic stats
    significant_hits = h3_summary %>%
        filter(grepl("condition",variable)) %>%
        mutate(ccf1 = str_extract(variable, "^ccf1(\\w+)", group = 1)) %>%
        filter(fdr<0.05)
    h3_mean_prop = h3 %>%
        group_by(h3, ccf1) %>%
        summarize(mean_prop = mean(prop))    
    all_hits = h3_summary %>%
        filter(grepl("condition",variable)) %>%
        mutate(ccf1 = str_extract(variable, "^ccf1(\\w+)", group = 1)) %>%
        mutate(h3 = cell_type) %>%
        left_join(h3_mean_prop)
    write_csv(all_hits, "analysis/differential_composition.csv")
    significant_hits_02 = all_hits %>%
        filter(fdr<0.05 | (fdr<0.2 & mean_prop>1e-3))
    write_csv(significant_hits_02, "analysis/differential_composition_significant.csv")
    
    significant_hits_02 # 21 significant interactions, 94 at 0.2 threshold, 29 at 0.2 thr + av prop
    dplyr::count(significant_hits_02, cell_type) # 16 types, 46 at 0.2 threshold, 22 at 0.2 thr + av prop
    dplyr::count(significant_hits_02, ccf1) # 9 regions, 18 at 0.2 threshold, 9 at 0.2 thr + av prop
    
    # MA plot
    h3_mean_prop = h3 %>%
        group_by(h3, ccf1) %>%
        summarize(mean_prop = mean(prop))    
    h3_summary %>%
        filter(grepl("condition",variable)) %>%
        mutate(ccf1 = str_extract(variable, "^ccf1(\\w+)", group = 1)) %>%
        mutate(var = paste(cell_type, ccf1, sep = "|"), h3=cell_type) %>%
        left_join(h3_mean_prop) %>%
        mutate(fdr_bin = cut(sign(Estimate)*fdr, breaks = c(-1,-0.2,-0.05,0,0.05,0.2,1))) %>%
        ggplot(aes(x=mean_prop,y=Estimate, label=var, color = fdr_bin)) +
        geom_point(size=0.5, show.legend = FALSE) +
        geom_hline(yintercept = 0, linetype="dashed") + 
        ggrepel::geom_text_repel(data= . %>% filter(fdr<0.05 | (fdr<0.2 & mean_prop>1e-3)), size=5, show.legend = FALSE, max.overlaps = Inf) +
        theme_classic(base_size = 16) +
        labs(x = "Average proportion", y="log2(FC)", col="FDR") +
        scale_x_log10() +
        scale_color_manual(values = c("gray70","skyblue3","blue3","brown1","darkorange","gray70"))
    ggsave("fig/diff_h3_ma.pdf")
    h3_summary %>%
        filter(grepl("condition",variable)) %>%
        mutate(ccf1 = str_extract(variable, "^ccf1([^:]+)", group = 1)) %>%
#        mutate(ccf1 = fct_reorder(ccf1, log(fdr)*((fdr<0.05)+0.01*(fdr>=0.05)), .fun = sum)) %>%
        mutate(ccf1 = fct_reorder(ccf1, log(fdr)/fdr, .fun = mean)) %>%
        ggplot(aes(x=ccf1,y=sign(Estimate)*(-log10(fdr)))) +
        geom_hline(yintercept = c(log10(0.05),-log10(0.05)), linetype="dashed") +
        geom_boxplot() +
        theme_classic() +
        theme(axis.text.x = element_text(angle=90,hjust=1,vjust=0.5)) +
        labs(x=NULL, y="Signed significance")
    sig_cols = c("SIG_DOWN"="blue3","WSIG_DOWN"="skyblue3","WSIG_UP"="darkorange","SIG_UP"="brown1")
    h3_summary %>%
        filter(grepl("condition",variable)) %>%
        mutate(ccf1 = str_extract(variable, "^ccf1([^:]+)", group = 1)) %>%
        filter(fdr<0.2) %>%
        mutate(significance = paste0(ifelse(fdr<0.05, "SIG", "WSIG"), "_", ifelse(Estimate>0, "UP", "DOWN"))) %>%
        mutate(significance = factor(significance, rev(names(sig_cols)))) %>%
        mutate(ccf1 = fct_reorder(ccf1, significance, .fun = length, .desc = TRUE)) %>%
#        group_by(ccf1, significance) %>%
        ggplot(aes(x=ccf1,fill=significance)) +
        geom_bar(show.legend = FALSE) +
        theme_classic(base_size = 16) +
        theme(axis.text.x = element_text(angle=90,hjust=1,vjust=0.5)) +
        labs(x=NULL, y="# cell types affected") +
        scale_fill_manual(values = sig_cols)
    ggsave("fig/diff_h3_by_ccf.pdf")
    h3_summary %>%
        filter(grepl("condition",variable)) %>%
        mutate(ccf = str_extract(variable, "^ccf1([^:]+)", group = 1)) %>%
        mutate(h2 = h3_to_h2[cell_type]) %>%
#        mutate(h2 = fct_reorder(h2, log(fdr)*((fdr<0.05)+0.1*(fdr>=0.05)), .fun = sum)) %>%
        mutate(h2 = fct_reorder(h2, log(fdr)/fdr, .fun = mean)) %>%
        ggplot(aes(x=h2,y=sign(Estimate)*(-log10(fdr)))) +
        geom_hline(yintercept = c(log10(0.05),-log10(0.05)), linetype="dashed") +
        geom_boxplot() +
        theme_classic(base_size = 16) +
        theme(axis.text.x = element_text(angle=90,hjust=1,vjust=0.5)) +
        labs(x=NULL, y="Signed significance")
    ggsave("fig/diff_h3_by_h2.pdf")
    
    h3_significant = h3_summary %>%
        filter(grepl("condition",variable) & fdr<0.05) %>%
        mutate(ccf = str_extract(variable, "^ccf1([^:]+)", group = 1)) %>%
        mutate(var = paste(cell_type, ccf, sep = "|")) %>%
        mutate(var = fct_reorder(var, fdr))
    h3 %>%
        mutate(var = paste(h3, ccf1, sep="|")) %>%
        filter(var %in% filter(h3_significant, Estimate>0)$var) %>%
        mutate(var = factor(var, levels(h3_significant$var))) %>%
        ggplot(aes(x=litter,y=prop,fill=condition)) +
        geom_col(position = "dodge") +
        facet_wrap(~ var) +
        theme_classic(base_size=12) +
        labs(x = "Litter", y="Proportion of cortical excitatory neurons",title="Enriched in ENU")
    h3 %>%
        mutate(var = paste(h3, ccf1, sep="|")) %>%
        filter(var %in% filter(h3_significant, Estimate<0)$var) %>%
        mutate(var = factor(var, levels(h3_significant$var))) %>%
        ggplot(aes(x=litter,y=prop,fill=condition)) +
        geom_col(position = "dodge") +
        facet_wrap(~ var) +
        theme_classic(base_size=12) +
        labs(x = "Litter", y="Proportion of cortical excitatory neurons",title="Enriched in CTL")
}

aggregate_ct_ccf = function(my_data, ct_var, ccf_var, min_total=0, min_max=0, min_min=0) {
    my_data %>%
        dplyr::count(sample, .data[[ct_var]], .data[[ccf_var]]) %>%
        complete(sample,.data[[ct_var]],.data[[ccf_var]], fill = list(n=0)) %>%
        group_by(.data[[ccf_var]], .data[[ct_var]]) %>%
        filter(sum(n)>=min_total & max(n) >= min_max & min(n) >= min_min) %>%
        group_by(sample) %>%
        mutate(prop = n / sum(n)) %>%
        separate(sample, c("litter","condition"), "\\_") %>%
        mutate(condition = ifelse(condition %in% c("1L","3L"), "CTL","ENU"))
}

summarize_varexp = function(lm_results) {
    lapply(lm_results, function(h) { as_tibble(anova(h),rownames = "variable") }) %>%
        bind_rows(.id = "cell_type") %>%
        mutate(pval = `Pr(>F)`, fdr = p.adjust(pval, "fdr")) %>%
        arrange(pval)
}

summarize_results = function(lm_results) {
    lapply(lm_results, function(h) { as_tibble(summary(h)$coefficients,rownames = "variable") }) %>%
        bind_rows(.id = "cell_type") %>%
        mutate(pval = `Pr(>|t|)`, fdr = p.adjust(pval, "fdr")) %>%
        arrange(pval)
}

if (sys.nframe() == 0) {
    main()
}