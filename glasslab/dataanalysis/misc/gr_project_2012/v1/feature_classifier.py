'''
Created on Jul 4, 2012

@author: karmel

The goal here is to see if we can tease out which features are 
predictive of characteristics of interest, such as pausing ratio
or transrepression. Worth a stab.
'''
from __future__ import division
from glasslab.dataanalysis.misc.gr_project_2012.v1.width_buckets import get_data_with_bucket_score
from glasslab.dataanalysis.machinelearning.logistic_classifier import LogisticClassifier
from glasslab.dataanalysis.misc.gr_project_2012.v1.elongation import get_rep_string,\
    total_tags_per_run
import os
import sys
import math
from sklearn.metrics.metrics import confusion_matrix
from random import shuffle

if __name__ == '__main__':
    learner = LogisticClassifier()
    dirpath = 'karmel/Desktop/Projects/Classes/Rotations/Finland_2012/GR_Project/classification_low_basal'
    dirpath = learner.get_path(dirpath)
    
    '''
    # First time file setup
    data = get_data_with_bucket_score(learner, dirpath)
    # Get mean of all values to compress
    grouped_prelim = data.groupby('glass_transcript_id',as_index=False)
    grouped = grouped_prelim.mean()
    
    
    # Except for tag_count which should be summed.
    grouped['tag_count'] = grouped_prelim['tag_count'].sum()['tag_count']
    
    
    grouped.to_csv(learner.get_filename(dirpath,'feature_vectors.txt'),sep='\t',index=False)
    
    
    raise Exception
    grouped = learner.import_file(learner.get_filename(dirpath, '../transcript_vectors.txt'))
    grouped = grouped.drop(['chr_name','ucsc_link_nod','gene_names',
                            'transcription_start','transcription_end'],axis=1)
    
    '''
    all_grouped = learner.import_file(learner.get_filename(dirpath, 'feature_vectors.txt'))
    
    # Filter so that we're only looking at genes with low basal expression.
    
    if False:
        # Can we predict pausing ratio?
        
        # Minimal ratio in KLA+Dex vs. KLA pausing
        try: min_ratio= float(sys.argv[1])
        except IndexError: min_ratio = 3
        try: force_choice = sys.argv[2].lower() == 'force'
        except IndexError: force_choice = False
        try: extra_dir = sys.argv[3]
        except IndexError: extra_dir = ''
        
        subdir = learner.get_filename(dirpath, 
                    extra_dir, 'pausing_ratio_{0}'.format(str(min_ratio).replace('.','_')))
        if force_choice: subdir = subdir + '_forced_choice'
        
        if not os.path.exists(subdir): os.makedirs(subdir)
        
        totals = total_tags_per_run()   
        for replicate_id in ('', 1, 2, 3, 4):
            rep_str = get_rep_string(replicate_id)
            
            grouped = all_grouped
            
            pausing_states = grouped.filter(regex=r'(kla_dex_\d_bucket_score|kla_dex_bucket_score|' +\
                                            r'kla_\d_bucket_score|kla_bucket_score)').fillna(0)
            
            # First do some column cleanup..
            # Remove cols that are the kla_dex_bucket target; or that are not from the current replicate;
            # or that are id fields
            regex = r'(?!(kla_dex_\d_bucket_score|kla_dex_bucket_score|.*gene.*tags))' +\
                            r'(({1}pu_1_kla)|{0})'.format('(?!(.*_\d_.*|.*_id))', 
                                                          rep_str and r'.*' + rep_str + r'.*|' or '')
            dataset = grouped.filter(regex=regex)
        
            
            labels = (pausing_states['kla_dex_{0}bucket_score'.format(rep_str)] > 0) & \
                        (pausing_states['kla_dex_{0}bucket_score'.format(rep_str)] \
                                /pausing_states['kla_{0}bucket_score'.format(rep_str)] >= min_ratio)
            
            '''
            #For checking non-trivials,
            # where Dex+KLA start tags are not equal to KLA start tags
            # (since that would make the kla_dex_gene_body_lfc
            # perfectly predictive.
            mask = (grouped['kla_dex_{0}gene_start_tags'.format(rep_str)]\
                        /grouped['kla_{0}gene_start_tags'.format(rep_str)])
            mask = mask.fillna(0)
            mask = mask.apply(lambda x: bool(x and abs(math.log(x,2)) > .5))
            
            labels = map(lambda x: int(x[0] and x[1]), zip(labels, mask))
            '''
            # Nulls should be zeroes. Fill those first.
            dataset = dataset.fillna(0)
            dataset = learner.normalize_data(dataset)
            
            if force_choice:
                # Skip if not force_choice because it takes too long on my laptop.
                
                # Now check reliability on non-trivial examples,
                # where Dex+KLA start tags are not equal to KLA start tags
                # (since that would make the kla_dex_gene_body_lfc
                # perfectly predictive.
                
                norm_factor = totals['kla_dex'][replicate_id or 0]/totals['kla'][replicate_id or 0]
                kla_dex_starts = grouped['kla_dex_{0}gene_start_tags'.format(rep_str)]
                kla_starts = grouped['kla_{0}gene_start_tags'.format(rep_str)]
                kla_dex_starts = kla_dex_starts.apply(lambda x: max(x, 1))
                kla_starts = kla_starts.apply(lambda x: max(x, 1))
                start_ratio = kla_dex_starts/(kla_starts*norm_factor)
                mask = start_ratio.apply(lambda x: abs(math.log(x,2))) > .25 # Approx 20% up or down
                print sum(mask)
                
                # We want a random subset of the non-trivial cases as our tests
                test_indices = dataset.ix[mask].index.values
                shuffle(test_indices)
                test_indices = test_indices[:int(len(test_indices)/3)]
                # Now test_indices is a list of indices we want. 
                test_vectors = dataset.ix[test_indices]
                test_labels = labels.ix[test_indices]
                print sum(test_labels)
                
                dataset = dataset[~dataset.index.isin(test_indices)].reset_index()
                labels = labels[~dataset.index.isin(test_indices)].reset_index()[0]
                
            print 'Total positive examples: ', sum(labels)
            
            classifier_type='logistic'
            best_err = 1.0
            best_c = 0
            best_chosen = []
            possible_k = [20, 10, 5, 2]
            mods = []
            for k in possible_k:
                if force_choice:
                    chosen = ['kla_{0}gene_body_lfc'.format(rep_str), 'dex_over_kla_{0}gene_body_lfc'.format(rep_str)]
                    
                else:
                    chosen = learner.get_best_features(dataset, labels, k=k)
                    
                num_features = len(chosen)
                
                err, c, mod = learner.run_nested_cross_validation(dataset, labels, columns=chosen,
                            classifier_type=classifier_type,
                            draw_roc=True, draw_decision_boundaries=force_choice,
                            title_suffix=replicate_id and 'Group {0}'.format(replicate_id) or 'Overall',
                            save_path_prefix=learner.get_filename(subdir,'plot_{0}group'.format(rep_str)),
                        )
                
                if err < best_err:
                    best_err = err
                    best_c = c
                    best_chosen = chosen
                    best_mod = mod
                if force_choice: break
        
            print "Best number of features: ", len(best_chosen)
            print "Best features: ", best_chosen
            print "Best C, MSE: ", best_c, best_err
            
            if force_choice:
                #mod = learner.get_model(classifier_type=classifier_type, C=best_c)
                #fitted = mod.fit(training_data, training_labels)
                test_vectors = test_vectors[best_chosen]
                predicted_probs = best_mod.predict_proba(test_vectors)
                err = learner.mse(predicted_probs, test_labels.values)
                print err
                print confusion_matrix(test_labels.values, predicted_probs[:,1] > .5)
                 
                learner.draw_roc(label_sets=[(test_labels.values, predicted_probs)], 
                                 save_path=learner.get_filename(subdir,'check_nontrivial_{0}group'.format(rep_str)))
                
                
                learner.draw_decision_boundaries(best_mod, best_chosen, 
                  test_vectors.as_matrix(), 
                  test_labels.values,
                  title = 'Decision Boundaries: ' + (replicate_id and 'Group {0}'.format(replicate_id) or 'Overall'), 
                  force_lim = [-3,3,-3,3],
                  save_path = learner.get_filename(subdir,'plot_{0}group'.format(rep_str))\
                                         + '_check_non_trivial_decision_boundaries.png'
                  )
    if True:
        # How about transrepression or derepression?
        try: force_choice = sys.argv[1].lower() == 'force'
        except IndexError: force_choice = False
        try: extra_dir = sys.argv[2]
        except IndexError: extra_dir = ''
        
        subdir = learner.get_filename(dirpath, 
                    'derepressed' + '_' + extra_dir)
        if force_choice: subdir = subdir + '_forced_choice'
        
        if not os.path.exists(subdir): os.makedirs(subdir)
        
        totals = total_tags_per_run()
        for replicate_id in ('', 1, 2, 3, 4):
            rep_str = get_rep_string(replicate_id)
            
            grouped = all_grouped
            grouped['dmso_{0}tags'.format(rep_str)] = grouped['dmso_{0}gene_start_tags'.format(rep_str)] +\
                                                        grouped['dmso_{0}gene_body_tags'.format(rep_str)]
            grouped['dmso_{0}tags_rpkm'.format(rep_str)] = grouped['dmso_{0}tags'.format(rep_str)]*(10**3*10**6)/grouped['length']\
                                                                /totals['dmso'][replicate_id or 0]
            #grouped = grouped[grouped['dmso_{0}tags_rpkm'.format(rep_str)] < .33]
            
            grouped = grouped[grouped['kla_{0}lfc'.format(rep_str)] >= 1]
            print len(grouped)
            
            dataset = grouped
            dataset['rep_{0}pausing_ratio'.format(rep_str)] = dataset['kla_dex_{0}bucket_score'.format(rep_str)] \
                                                            /dataset['kla_{0}bucket_score'.format(rep_str)]
            targets = dataset.reset_index().filter(regex=r'.*_lfc')
            # Get rid of lfc columns    
            dataset = dataset.reset_index().filter(regex=r'(?!(.*_lfc|dex_over_kla.*|.*id|.*_tags))')
            # And then columns from other replicates
            dataset = dataset.reset_index().filter(regex=r'(({1}pu_1_kla)|{0})'.format('(?!(.*_\d_.*|.*_id))', 
                                                          rep_str and r'.*' + rep_str + r'.*|' or ''))
            print dataset.columns
            
            labels = targets['dex_over_kla_{0}lfc'.format(rep_str)] <= -.58
            labels_2 = targets['kla_{0}lfc'.format(rep_str)] >= 1
            
            labels = map(lambda x: int(x[0] and x[1]), zip(labels, labels_2))
            print 'Total examples, positive: ', len(dataset), sum(labels)
            # Nulls should be zeroes. Fill those first.
            dataset = dataset.fillna(0)
            dataset = learner.make_numeric(dataset)
            dataset = learner.normalize_data(dataset)
            
            
            best_err = 1.0
            best_c = 0
            best_chosen = []
            possible_k = [20, 10, 5, 2]
            for k in possible_k:
                if force_choice:
                    chosen = ['gr_kla_dex_tag_count', 'gr_dex_tag_count'.format(rep_str)]
                    #chosen = ['kla_{0}bucket_score'.format(rep_str), 'tag_count']
                else:
                    chosen = learner.get_best_features(dataset, labels, k=k)
                    
                num_features = len(chosen)
                
                err, c, mod = learner.run_nested_cross_validation(dataset, labels, columns=chosen,
                            draw_roc=True, draw_decision_boundaries=force_choice,
                            title_suffix=replicate_id and 'Group {0}'.format(replicate_id) or 'Overall',
                            save_path_prefix=learner.get_filename(subdir,'plot_{0}group'.format(rep_str)),
                        )
                
                if err < best_err:
                    best_err = err
                    best_c = c
                    best_chosen = chosen
                    
                if force_choice: break
        
            print "Best number of features: ", len(best_chosen)
            print "Best features: ", best_chosen
            print "Best C, MSE: ", best_c, best_err

            
            
        