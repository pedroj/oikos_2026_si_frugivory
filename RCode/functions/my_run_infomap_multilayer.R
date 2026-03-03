my_run_infomap_multilayer<-
    function (M, infomap_executable = "Infomap", flow_model = NULL, 
              silent = T, trials = 100, seed = NULL, relax = F, multilayer_relax_rate = 0.1, 
              multilayer_relax_limit = NULL, multilayer_relax_limit_up = NULL, 
              multilayer_relax_limit_down = NULL, temporal_network = F, 
              run_standalone = T, remove_auxilary_files = T, ...) 
    {
        if (!"inter" %in% names(M)) {
            if (any(M$extended_ids$layer_from != M$extended_ids$layer_to) && 
                relax == FALSE) {
                M$inter <- as.data.frame(M$extended_ids[M$extended_ids$layer_from != 
                                                            M$extended_ids$layer_to, ])
                M$inter <- as.tbl(M$inter)
                M$inter <- M$inter %>% dplyr::mutate_all(as.numeric)
            }
        }
        else {
            inter <- NULL
        }
        intra <- NULL
        if (!"intra" %in% names(M) && relax == T) {
            intra <- M$extended_ids[M$extended_ids$layer_from == 
                                        M$extended_ids$layer_to, c("layer_from", "node_from", 
                                                                   "node_to", "weight")]
            colnames(intra)[1] <- "layer"
        }
        else if (!"intra" %in% names(M) && relax == F) {
            intra <- M$extended_ids[M$extended_ids$layer_from == 
                                        M$extended_ids$layer_to & as.numeric(M$extended_ids$weight) != 
                                        0, ]
        }
        if (!is.null(intra)) {
            M$intra <- as.data.frame(intra)
            M$intra <- M$intra %>% as.tbl(M$intra) %>% dplyr::mutate_all(as.numeric)
        }
        if (check_infomap(infomap_executable) == F) {
            stop("Error in Infomap stand-alone file.")
        }
        if (class(M) != "multilayer") {
            stop("M must be of class multilayer")
        }
        arguments <- paste("--tree -2 -N ", trials, sep = "")
        arguments <- ifelse(!is.null(seed), paste(arguments, "--seed", 
                                                  seed), arguments)
        arguments <- ifelse(!is.null(flow_model), paste(arguments, 
                                                        "-f", flow_model), arguments)
        arguments <- ifelse(silent, paste(arguments, "--silent"), 
                            arguments)
        arguments <- paste(arguments, ...)
        if (relax == F) {
            print("Using interlayer edge values to determine flow between layers.")
            write_lines("*Multilayer", "infomap_multilayer.txt")
            write_delim(M$intra, "infomap_multilayer.txt", delim = " ", 
                        append = T)
            write_delim(M$inter, "infomap_multilayer.txt", delim = " ", 
                        append = T)
        }
        else {
            if (ncol(M$intra) == 5) {
                stop("Cannot use relax rates with extended format of intralayer edges. See function create_multilayer_object.")
            }
            print("Using global relax to determine flow between layers.")
            write_lines("*Intra", "infomap_multilayer.txt")
            write_delim(M$intra, "infomap_multilayer.txt", delim = " ", 
                        append = T)
            if (!is.null(M$inter)) {
                if (ncol(M$inter) == 5) {
                    stop("Cannot use relax rates with extended format of interlayer edges. See function create_multilayer_object.")
                }
                print("Global relax will be constrained by interlayer edges.")
                write_lines("*Inter", "infomap_multilayer.txt", append = T)
                write_delim(M$inter, "infomap_multilayer.txt", delim = " ", 
                            append = T)
            }
            arguments <- ifelse(!is.null(multilayer_relax_rate), 
                                paste(arguments, "--multilayer-relax-rate", multilayer_relax_rate), 
                                arguments)
            arguments <- ifelse(!is.null(multilayer_relax_limit), 
                                paste(arguments, "--multilayer-relax-limit", multilayer_relax_limit), 
                                arguments)
            arguments <- ifelse(!is.null(multilayer_relax_limit_up), 
                                paste(arguments, "--multilayer-relax-limit-up", multilayer_relax_limit_up), 
                                arguments)
            arguments <- ifelse(!is.null(multilayer_relax_limit_down), 
                                paste(arguments, "--multilayer-relax-limit-down", 
                                      multilayer_relax_limit_down), arguments)
        }
        call <- paste("./", infomap_executable, " infomap_multilayer.txt . ", 
                      arguments, sep = "")
        if (run_standalone == T) {
            print(call)
            system(call)
        }
        else {
            print("Please run Infomap online at https://www.mapequation.org/infomap/ using the following arguments (copy-paste):")
            print(arguments)
            invisible(readline(prompt = "After running, download statenodes results and press [ENTER] when done"))
            if (!file.exists("network_states.tree")) {
                stop("Result file network_states.tree was not found. Did you download results?")
            }
            file.rename(from = "network_states.tree", to = "infomap_multilayer_states.tree")
        }
        L_output <- parse_number(read_lines("infomap_multilayer_states.tree")[6])
        modules <- suppressMessages(read_delim("infomap_multilayer_states.tree", 
                                               delim = " ", skip = 11, col_names = c("path", "flow", 
                                                                                     "name", "state_id", "node_id", "layer_id")))
        modules %<>% dplyr::filter(flow > 0) %>% dplyr::select(path, node_id, layer_id, 
                                                 flow) %>% separate(path, into = c("module", "leaf_id"), 
                                                                    sep = ":") %>% dplyr::mutate_at(.vars = 1:4, as.integer) %>% 
            dplyr::full_join(M$nodes, "node_id") %>% dplyr::select(node_id, starts_with("module"), 
                                                     everything(), -leaf_id) %>% dplyr::arrange(node_id, layer_id)
        if (any(is.na(modules$module))) {
            print("No modules were assigned nodes that lacked flow: ")
            print(modules[is.na(modules$module), c("node_id", "node_name")])
            modules <- modules[!is.na(modules$module), ]
        }
        if (temporal_network) {
            print("Reorganizing modules...")
            renamed_moduels <- modules %>% dplyr::distinct(module, layer_id) %>% 
                dplyr::arrange(module, layer_id)
            x <- c(1, table(renamed_moduels$module))
            module_birth_layers <- renamed_moduels %>% dplyr::slice(cumsum(x)) %>% 
                dplyr::arrange(layer_id, module)
            module_renaming <- data.frame(module = module_birth_layers$module, 
                                          module_renamed = 1:max(module_birth_layers$module))
            modules %<>% dplyr::left_join(module_renaming, "module") %>% 
                dplyr::select(-module) %>% dplyr::rename(module = module_renamed)
        }
        if (remove_auxilary_files) {
            print("Removing auxilary files...")
            file.remove("infomap_multilayer_states.tree")
            file.remove("infomap_multilayer.txt")
            file.remove("infomap_multilayer.tree")
        }
        print(paste("Partitioned into ", max(modules$module), " modules.", 
                    sep = ""))
        out <- list(call = call, L = L_output, m = max(modules$module), 
                    modules = modules)
        class(out) <- "infomap_multilayer"
        return(out)
    }