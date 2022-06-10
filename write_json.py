import json
import logging
import os
import __main__ as main

def write_json(data, file=None, outdir=None):
	if not file:
		if not outdir:
			logging.info(f'file or outdir should be defined !')
			return None
		else:
			file = os.path.join(outdir, os.path.splitext(os.path.split(os.path.realpath(main.__file__))[1])[0] +'.json')
	if type(data) == dict:
		if os.path.isfile(file):
			logging.info(f'write json: {file} exists')
			return None
		try:
			with open(file, 'w') as H:
				json.dump(data, H, indent=2)
			logging.info(f'write json to {file}')
			return True
		except Exception as e:
			logging.error(e)
			return None
	else:
		logging.info(f'write json: data is not dict type')
		return None




