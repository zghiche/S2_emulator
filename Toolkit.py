import yaml
import cppyy

with open('config.yaml', "r") as afile:
    params = yaml.safe_load(afile)

def define_map(params):
    enum = cppyy.gbl.l1thgcfirmware.Step
    map_custom = cppyy.gbl.std.map[enum, "unsigned int"]()

    for d in params['stepLatency']:
      for latency in d.keys():
        map_custom[int(vars(enum)[latency])] = d[latency]

    params['stepLatency'] = map_custom
