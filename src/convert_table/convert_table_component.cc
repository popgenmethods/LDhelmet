//this component was written by Jeffrey P. Spence in 2016
//adapted from table_gen written by Andrew H. Chan
// email: spence.jeffrey@berkeley.edu

#include "convert_table/convert_table_component.h"
#include "convert_table/convert_table_options.h"
#include <stddef.h>
#include <iostream>
#include <stdint.h>
#include <stdio.h>
#include <fstream>
#include <algorithm>
#include <cstdlib>
#include <unordered_map>
#include <vector>
#include "common/read_confs.h"
#include "common/rho_finder.h"
#include "common/version_number.h"
#include "table_gen/output_writer.h"
#include "table_gen/table_management.h"

std::vector<std::string> split(std::string const &input) {
    std::istringstream buffer(input);
    std::vector<std::string> ret;
    
    std::copy(std::istream_iterator<std::string>(buffer),
              std::istream_iterator<std::string>(),
              std::back_inserter(ret));
    return ret;
}


int ConvertTableComponent(std::string const base_command, int argc, char **argv) {
    //uses the same format as table_gen
    uint64_t version_number = TABLE_GEN_VERSION_OUTPUT;
    uint64_t version_bit_string = TABLE_GEN_SALT + version_number;
    std::string version_string = TABLE_GEN_VERSION_STRING;
    CmdLineOptionsConvertTable cmd_line_options(base_command,
                                            argc,
                                            argv,
                                            version_string);
    if (!cmd_line_options.success()) {
        std::exit(1);
    }
    printf("%s\n", ShowCommandLineOptions(argc, argv).c_str());
    std::string const &input_file = cmd_line_options.input_file_;
    std::string const &output_file = cmd_line_options.output_file_;
    std::string const &conf_file = cmd_line_options.config_file_;
    
    
    
    
    std::ifstream ldhatfile;
    ldhatfile.open(input_file);
    
    std::string curr_line;
    //read in header
    std::getline(ldhatfile, curr_line, ' '); //num haps
    size_t num_haps = std::stoi(curr_line);
    std::getline(ldhatfile, curr_line); //num configs
    logfactorials.push_back(0.0);
    for(size_t n = 1; n <= num_haps; n++){
        logfactorials.push_back(logfactorials[n-1] + log(n));
    }
    std::getline(ldhatfile, curr_line, ' '); //number that does nothing
    std::getline(ldhatfile, curr_line); //theta
    std::vector<double> thetas;
    thetas.push_back(std::stod(curr_line));
    
    //get rhos
    std::getline(ldhatfile, curr_line);
    std::vector<std::string> rho_array = split(curr_line);
    std::vector<double> rhoArg;
    if (rho_array.size() == 2){ //this is an LDHat table
        std::cout << "LDhat table\n" << std::endl;
        rhoArg.push_back(0.0);  //LDHat tables always start at 0.0
        rhoArg.push_back(std::stod(rho_array[1]) / (std::stod(rho_array[0]) - 1));   //stepsize is max_rho / num_rhos
        rhoArg.push_back(std::stod(rho_array[1]));
    } else {                    //this is LDHelmet format table
        for(int i = 0; i < rho_array.size(); i++){
            rhoArg.push_back(std::stod(rho_array[i]));
        }
    }
    boost::tuple<std::vector<double>, std::vector<double> > rho_grid = ParseRhoRange(rhoArg);
    std::vector<double> rhos = GetRhoList(rho_grid);
    assert(rhos.size() > 0);
    if (rhos.front() != 0.0) {
        fprintf(stderr,
                "The rho values must begin at 0.0.\n");
        std::exit(1);
    }
    
    
    
    //Print out what we've read in, so far
    printf("theta: ");
    for (int i = 0; i < static_cast<int>(thetas.size()); ++i) {
        printf("%g ", thetas[i]);
    }
    printf("\n");
    
    printf("rho values: ");
    for (int i = 0; i < static_cast<int>(rhos.size()); ++i) {
        printf("%g ", rhos[i]);
    }
    printf("\n");
    
    std::set<size_t> degrees;
    degrees.insert(2 * num_haps);
    
    //Check if there's a config file, and if it has missing data
    std::vector<Conf> missing_confs;
    if(conf_file.size() > 0){
        std::cout << "Configuration File = " << conf_file << std::endl;
        std::vector<Conf> original_conf_list = LoadConfigurations(conf_file);
        for(size_t conf_idx = 0; conf_idx < original_conf_list.size(); conf_idx++){
            if(original_conf_list[conf_idx].ComputeM() < num_haps){
                missing_confs.push_back(original_conf_list[conf_idx]);
                degrees.insert(original_conf_list[conf_idx].ComputeDegree());
            }
            if(original_conf_list[conf_idx].ComputeSize() > num_haps){
                fprintf(stderr,
                        "There is a config in your sample with more samples than in your table.\n");
                std::exit(1);
            }
        }
    }
    
    //Read in all of the configs
    printf("Reading likelihoods from file.\n");
    std::vector<std::pair<Conf, std::vector<double> > > conf_list_master;
    while(std::getline(ldhatfile,curr_line)){
        if(curr_line.compare("") == 0){
            continue;
        }
        std::vector<std::string> line_array = split(curr_line);
        //read in the number of the different types
        int aa = std::stoi(line_array[2]);
        int ab = std::stoi(line_array[3]);
        int ba = std::stoi(line_array[4]);
        int bb = std::stoi(line_array[5]);
        
        int rho_offset = 0;
        while (line_array[rho_offset].compare(":") != 0 && rho_offset < line_array.size()){
            rho_offset += 1;
        }
        if (rho_offset == line_array.size()){
            fprintf(stderr, "Config likelihood lines must contain ':' \n");
            std::exit(1);
        }
        if (rho_offset + rhos.size() + 1 != line_array.size()){
            fprintf(stderr, "Number of table entries does not match number of rhos \n");
            std::exit(1);
        }
        std::vector<double> rhoLogLiks;
        for(int rhoIdx = rho_offset + 1; rhoIdx < line_array.size(); ++rhoIdx)
        {
            double currLoglik = exp(std::stod(line_array[rhoIdx]));
            rhoLogLiks.push_back(currLoglik);
        }
        
        //undo the symmetries that LDHat assumes and LDHelmet does not
        Conf currConf1(0, 0, 0, 0, aa, ab, ba, bb);
        std::pair<Conf, std::vector<double> > toAdd1(currConf1, rhoLogLiks);
        conf_list_master.push_back(toAdd1);
        
        //swap allele at first locus
        Conf currConf2(0, 0, 0, 0, ba, bb, aa, ab);
        std::pair<Conf, std::vector<double> > toAdd2(currConf2, rhoLogLiks);
        conf_list_master.push_back(toAdd2);
        
        //swap allele at second locus
        Conf currConf3(0, 0, 0, 0, ab, aa, bb, ba);
        std::pair<Conf, std::vector<double> > toAdd3(currConf3, rhoLogLiks);
        conf_list_master.push_back(toAdd3);
        
        //swap allele at both loci
        Conf currConf4(0, 0, 0, 0, bb, ba, ab, aa);
        std::pair<Conf, std::vector<double> > toAdd4(currConf4, rhoLogLiks);
        conf_list_master.push_back(toAdd4);
        
        //swap loci
        Conf currConf5(0, 0, 0, 0, aa, ba, ab, bb);
        std::pair<Conf, std::vector<double> > toAdd5(currConf5, rhoLogLiks);
        conf_list_master.push_back(toAdd5);
        
        //swap loci, swap allele at first locus
        Conf currConf6(0, 0, 0, 0, ba, aa, bb, ab);
        std::pair<Conf, std::vector<double> > toAdd6(currConf6, rhoLogLiks);
        conf_list_master.push_back(toAdd6);
        
        //swap loci, swap allele at second locus
        Conf currConf7(0, 0, 0, 0, ab, bb, aa, ba);
        std::pair<Conf, std::vector<double> > toAdd7(currConf7, rhoLogLiks);
        conf_list_master.push_back(toAdd7);
        
        //swap loci, both alleles
        Conf currConf8(0, 0, 0, 0, bb, ab, ba, aa);
        std::pair<Conf, std::vector<double> > toAdd8(currConf8, rhoLogLiks);
        conf_list_master.push_back(toAdd8);
        
    }
    
    
    
    if(conf_file.size() > 0){
        printf("Marginalizing to obtain configs with missing data (if any).  This may take a while.\n");
        
        
        //first fill the hash_map with the fully specified configs.
        std::unordered_map<std::string, std::vector<double>> like_map;
        for(size_t conf_idx = 0; conf_idx < conf_list_master.size(); conf_idx++){
            like_map[conf_list_master[conf_idx].first.GetString()] = conf_list_master[conf_idx].second;
        }
        //now compute the missing configs.
        for(size_t miss_idx = 0; miss_idx < missing_confs.size(); miss_idx++){
            Conf this_conf = missing_confs[miss_idx];
            std::vector<double> these_likelihoods = marginalProbs(this_conf, like_map);
            //check that these are good
            for(size_t i = 0; i < these_likelihoods.size(); i++){
                if(these_likelihoods[i] > 1.0 || these_likelihoods[i] < 0){
                    std::cout << "Encountered a wrong probability:" << std::endl;
                    std::cout << log(these_likelihoods[0]) << std::endl;
                }
            }
            std::pair<Conf, std::vector<double> > toAdd(this_conf, these_likelihoods);
            conf_list_master.push_back(toAdd);
        }
                
    }
    
    struct sort_pre {
        bool operator()(const std::pair<Conf,std::vector<double > > &left, const std::pair<Conf,std::vector<double > > &right) {
            return left.first < right.first;
        }
    };
    std::sort(conf_list_master.begin(), conf_list_master.end(), sort_pre());
    std::vector<Conf> conf_list;
    for(int i = 0; i < conf_list_master.size(); ++i){
        conf_list.push_back(conf_list_master[i].first);
    }
    std::vector<size_t> degree_seps;
    uint32_t max_degree, max_sample_size, max_locus;
    boost::tie(degree_seps, max_degree, max_sample_size, max_locus) =
                    PreprocessConfs(conf_list);
    
    
    
    //Start writing:
    InputConfBinaryWriter input_conf_binary_writer(output_file,
                                                   conf_list,
                                                   degree_seps);
    {
        // Write input confs.
        InputConfBinaryWriter input_conf_binary_writer(output_file,
                                                       conf_list,
                                                       degree_seps);
        {
            // Write version number.
            int num_written;
            num_written = fwrite(reinterpret_cast<void const *>(&version_bit_string),
                                 sizeof(version_bit_string),
                                 1, input_conf_binary_writer.fp_);
            assert(num_written == 1);
            
            // Write num confs.
            uint64_t num_conf_list = conf_list.size();
            num_written = fwrite(reinterpret_cast<void const *>(&num_conf_list),
                                 sizeof(num_conf_list),
                                 1, input_conf_binary_writer.fp_);
            assert(num_written == 1);
            
            // Write theta.
            double theta = thetas.front();
            num_written = fwrite(reinterpret_cast<void const *>(&theta),
                                 sizeof(theta),
                                 1, input_conf_binary_writer.fp_);
            assert(num_written == 1);

            // This turns on interpolation
            uint8_t interpolate = 1;
            num_written = fwrite(reinterpret_cast<void const *>(&interpolate),
                                 sizeof(interpolate),
                                 1, input_conf_binary_writer.fp_);
            assert(num_written == 1);

            
            // Write number of rho segments, excluding end point.
            if (rho_grid.get<0>().size() == 0) {
                fprintf(stderr,
                        "Error: The size of rho_grid is 0.\n");
                std::exit(1);
            }
            uint64_t num_rho_segments = rho_grid.get<0>().size() - 1;
            num_written = fwrite(reinterpret_cast<void const *>(&num_rho_segments),
                                 sizeof(num_rho_segments),
                                 1, input_conf_binary_writer.fp_);
            assert(num_written == 1);
            
            // Write rho segments.
            assert(rho_grid.get<0>().size() > 0);
            for (uint32_t i = 0; i < rho_grid.get<0>().size() - 1; ++i) {
                double start = rho_grid.get<0>()[i];
                double delta = rho_grid.get<1>()[i];
                num_written = fwrite(reinterpret_cast<void const *>(&start),
                                     sizeof(start),
                                     1, input_conf_binary_writer.fp_);
                assert(num_written == 1);
                
                num_written = fwrite(reinterpret_cast<void const *>(&delta),
                                     sizeof(delta),
                                     1, input_conf_binary_writer.fp_);
                assert(num_written == 1);
            }
            num_written = fwrite(reinterpret_cast<void const *>(
                                                                &(rho_grid.get<0>().back())),
                                 sizeof(rho_grid.get<0>().back()),
                                 1, input_conf_binary_writer.fp_);
            assert(num_written == 1);
            
            // Write confs.
            for (std::vector<Conf>::const_iterator citer = conf_list.begin();
                 citer != conf_list.end();
                 ++citer) {
                Conf::BinaryRep::RepType rep = Conf::BinaryRep(*citer).rep();
                num_written = fwrite(reinterpret_cast<void const *>(&rep),
                                     sizeof(rep),
                                     1, input_conf_binary_writer.fp_);
                assert(num_written == 1);
            }
        }
        
        assert(thetas.size() == 1);
        printf("Writing table\n");
        for (size_t theta_id = 0; theta_id < thetas.size(); ++theta_id) {
            for (size_t rho_id = 0; rho_id < rhos.size(); ++rho_id) {
                std::cout << "Working on Rho ID:" << rho_id << std::endl;
                
                Vec8 table(max_degree + 1);
                std::set<size_t>::iterator it;
                for (it = degrees.begin(); it != degrees.end(); it++) {
                    AllocateMemoryToTable(&table, *it);
                    for (int confIdx = 0; confIdx < conf_list_master.size(); ++confIdx){
                        if(conf_list_master[confIdx].first.ComputeDegree() == *it){
                            SetTable(&table, conf_list_master[confIdx].first, conf_list_master[confIdx].second[rho_id]);
                        }
                    }
                    input_conf_binary_writer.Write(*it, table);
                }
            }
        }
    }
    return 0;
}


std::vector<double> marginalProbs(Conf config, std::unordered_map<std::string, std::vector<double>> &like_map){
    std::string config_string = config.GetString();
    
    if(like_map.count(config_string) == 1){
        return like_map[config_string];
    }
    //this config is not in the map.  First we move toward being fully specified, and then we move toward having more configs
    std::vector<double> to_return;
    Conf temp_config;
    if(config.Geta(0) > 0){
        //remove a 0_
        temp_config = config.Adda(0, -1);
        
        //add a 00
        temp_config = temp_config.Addr(0, 0, 1);
        std::vector<double> prob1 = marginalProbs(temp_config, like_map);
        temp_config = temp_config.Addr(0, 0, -1);
        
        //add a 01
        temp_config = temp_config.Addr(0, 1, 1);
        std::vector<double> prob2 = marginalProbs(temp_config, like_map);
        temp_config = temp_config.Addr(0, 1, -1);
        
        temp_config = temp_config.Adda(0, 1);
        to_return = vectorSum(prob1, prob2);
        like_map[config_string] = to_return;
        return marginalProbs(config, like_map);
    }
    if(config.Geta(1) > 0){
        //remove a 1_
        temp_config = config.Adda(1, -1);
        
        //add a 10
        temp_config = temp_config.Addr(1, 0, 1);
        std::vector<double> prob1 = marginalProbs(temp_config, like_map);
        temp_config = temp_config.Addr(1, 0, -1);

        
        //add a 11
        temp_config = temp_config.Addr(1, 1, 1);
        std::vector<double> prob2 = marginalProbs(temp_config, like_map);
        temp_config = temp_config.Addr(1, 1, -1);
        
        temp_config = temp_config.Adda(1, 1);
        to_return = vectorSum(prob1, prob2);
        like_map[config_string] = to_return;
        return marginalProbs(config, like_map);
    }
    if(config.Getb(0) > 0){
        
        //remove a _0
        temp_config = config.Addb(0, -1);
        
        //add a 00
        temp_config = temp_config.Addr(0, 0, 1);
        std::vector<double> prob1 = marginalProbs(temp_config, like_map);
        temp_config = temp_config.Addr(0, 0, -1);
        
        //add a 10
        temp_config = temp_config.Addr(1, 0, 1);
        std::vector<double> prob2 = marginalProbs(temp_config, like_map);
        temp_config = temp_config.Addr(1, 0, -1);
        
        temp_config = temp_config.Addb(0, 1);
        to_return = vectorSum(prob1, prob2);
        like_map[config_string] = to_return;
        return marginalProbs(config, like_map);
    }
    if(config.Getb(1) > 0){
        
        //remove a _1
        temp_config = config.Addb(1, -1);

        //add a 01
        temp_config = temp_config.Addr(0, 1, 1);
        std::vector<double> prob1 = marginalProbs(temp_config, like_map);
        temp_config = temp_config.Addr(0, 1, -1);

        //add a 11
        temp_config = temp_config.Addr(1, 1, 1);
        std::vector<double> prob2 = marginalProbs(temp_config, like_map);
        temp_config = temp_config.Addr(1, 1, -1);
        
        temp_config = temp_config.Addb(1, 1);
        to_return = vectorSum(prob1, prob2);
        
        like_map[config_string] = to_return;
        return marginalProbs(config, like_map);
    }
    
    //If we are here, the config is fully specified, so we compute it by marginalization over k+1st ind
    
    for(size_t a_allele = 0; a_allele < 2; a_allele++){
        for(size_t b_allele = 0; b_allele < 2; b_allele++){
            temp_config = config.Addr(a_allele, b_allele, 1);
            to_return = vectorSum(to_return, marginalProbs(temp_config, like_map));
        }
    }
    like_map[config_string] = to_return;
    return marginalProbs(config, like_map);
}

// returns x+y for vectors x and y.  if either x or y is empty, it's assumed to be the zero vector
std::vector<double> vectorSum(std::vector<double> vec1, std::vector<double> vec2){
    if(vec1.size() == 0){
        assert (vec2.size() > 0);
        return vec2;
    }
    if(vec2.size() == 0){
        assert (vec1.size() > 0);
        return vec1;
    }
    assert (vec1.size() == vec2.size());
    std::vector<double> to_return;
    for (size_t i = 0; i < vec1.size(); i++){
        to_return.push_back(vec1[i] + vec2[i]);
    }
    return to_return;
}
