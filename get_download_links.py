#!/usr/bin/env python

import requests, sys, re
base = "https://zenodo.org"
args = sys.argv
if len(args)==2:
	if args[1][:len(base)] != base:
		print("Invalid argument, please provide a zenodo record URL")
	else:
		response = requests.get(args[1])

		filename_pat = re.compile(r'(?<=<a class="filename" href=")(.*?)(?=">)(?:">)(?P<filename>(?<=">).*?(?=</a>))')
		matches = filename_pat.findall(str(response.content))

		links, names = zip(*matches)
		links = map(lambda x: base+x, list(links))
		
		for l, n in zip(links, names):
			print(l,n)

