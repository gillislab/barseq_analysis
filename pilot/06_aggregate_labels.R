
library(tidyverse)


main = function() {
    class_labels = read_labels("whole") %>%
        rename(class=cluster_name) %>%
        mutate(class = simplify_class_name(class))
    glu_labels = read_labels("subglu")
    gaba_labels = read_labels("gaba")
    
    labels = left_join(class_labels, bind_rows(glu_labels, gaba_labels))
    write_csv(labels, "analysis/labels.csv")
}

read_labels = function(subdir) {
    input_dir = file.path("analysis", subdir)
    clusters = read_csv(file.path(input_dir, "cluster.csv"))
    annotation = read_csv(file.path(input_dir, "cluster_annotation.csv"))
    result = inner_join(clusters, annotation, by = c(label = "cluster_id")) %>%
        select(-label)
    return(result)
}

simplify_class_name = function(label) {
    result = label
    result[startsWith(label, "GLU")] = "Glutamatergic"
    result[startsWith(label, "GABA")] = "GABAergic"
    return(result)
}

if (sys.nframe() == 0) {
    main()
}
