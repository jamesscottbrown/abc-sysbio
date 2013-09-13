import numpy, multiprocessing

import cudasim.Lsoda_mg as Lsoda
import cudasim.EulerMaruyama_mg as EulerMaruyama
import cudasim.Gillespie_mg as Gillespie

class cuda_model:
    
    # instantiation
    def __init__(self, name, nspecies, nparameters, prior, x0prior, source, integration, fit, dt, beta, timepoints, logp, ngpu):
        self.nspecies = nspecies
        self.name = name

        ## combine the parameters with the species
        self.kparameters = nparameters
        self.nparameters = nparameters + nspecies
        
        self.prior = [x[:] for x in prior]
        for x in x0prior:
            self.prior.append( x[:] )
        
        self.source = source
        self.integration = integration
        self.fit = fit
        self.cudaCode = self.name +  '.cu' 
        self.dt = dt
        self.beta = beta
        self.timepoints = timepoints
        self.logp = logp
        self.ngpu = ngpu
       
    def simulate(self, p, t, n, beta):
        gpu_threads = []
        output_cpu  = multiprocessing.Queue()

        # distribute the particles across the cards
        n_per_card = [0 for i in range(self.ngpu)]
        nc = int(n/self.ngpu)
        for i in range(self.ngpu-1):
            n_per_card[i] = int(n/self.ngpu)

        n_per_card[self.ngpu-1] = n - (self.ngpu-1)*nc

        ## print n, n_per_card, numpy.shape(p)
        ## card_ids = [0, 1, 3, 4]
   
        for c in range(self.ngpu):

            ### Create local parameter and species arrays
            species = numpy.zeros([n_per_card[c],self.nspecies])
            pp = numpy.zeros([n_per_card[c],self.kparameters])

            ### Fill species and parameters
            for i in range(n_per_card[c]):
                place_mark = sum(n_per_card[0:c])
                species[i,:] = p[place_mark + i][self.kparameters:self.nparameters]
            
                if self.logp == False:
                    pp[i,:] = p[place_mark + i][0:self.kparameters]
                else:
                    pp[i,:] = numpy.power(10,p[place_mark + i][0:self.kparameters])

            ### Run on multiple GPUs
            gpu_thread = Lsoda_mg.Lsoda(self.timepoints, self.cudaCode, pp, species, output_cpu, card=c, dt=self.dt, dump=False, info=False, timing=False)
            gpu_threads.append(gpu_thread)
            gpu_thread.start()

        result_dict = {}
        for gpu_pro in gpu_threads:
            id, results = output_cpu.get(gpu_pro)
            ## print id, np.shape(results)
            result_dict[id] = results

        ## print result_dict.keys()
        result = numpy.vstack( result_dict.values() )

        ## print numpy.shape(result)

        return result

