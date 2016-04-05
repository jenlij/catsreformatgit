//
//  Hash.h
//  CatsReformat
//
//  Created by Jennifer Johnson on 5/2/15.
//  Copyright (c) 2015 Jennifer Johnson. All rights reserved.
//

#ifndef __HASH_H__
#define __HASH_H__

unsigned int hash ( const unsigned char * key, unsigned int length, unsigned int initval);

unsigned int hash_no_case ( const unsigned char * key, unsigned int length, unsigned int initval);

#endif
