// ECE4122/6122 RSA Encryption/Decryption assignment
// Fall Semester 2015

#include <iostream>
#include "RSA_Algorithm.h"

using namespace std;

int main()
{
  // Instantiate the one and only RSA_Algorithm object
  RSA_Algorithm RSA;
  
  // Loop from sz = 32 to 1024 inclusive
  // for each size choose 10 different key pairs
  // For each key pair choose 10 differnt plaintext 
  // messages making sure it is smaller than n.
  // If not smaller then n then choose another
  // For eacm message encrypt it using the public key (n,e).
  // After encryption, decrypt the ciphertext using the private
  // key (n,d) and verify it matches the original message.

  // your code here
  mpz_class M; // Original message
  mpz_class C; // Encrypter M
  mpz_class M_result; // Decrypted C

  size_t iter_times = 100; // Set iteration times, 100 or 10?
  size_t max_sz = 1024;

  size_t sz;

  size_t true_count = 0;
  size_t false_count = 0;

  for (sz = 32; sz <= max_sz; sz *= 2){ 
  // sz is doubled every time
    for (size_t i = 0; i < iter_times; ++i){
      // Create 100 different public/private key pairs
      RSA.GenerateRandomKeyPair(sz);
      for (size_t j = 0; j < iter_times; ++j){
        // For each pair of key, create 100 random messages
        // Set size of message as sizeof(n)-1;
        
        do{
          M = (RSA.rng).get_z_bits(2 * sz - 1);
        }while(M >= RSA.n);
        
        //RSA.PrintNDE();
        // Randon message

        RSA.PrintM(M);
        C = RSA.Encrypt(M);
        M_result = RSA.Decrypt(C);

        
       // cout << sz << "-" << i << "-" << j << endl;

        // if (M == M_result) {
        //   cout << "true" << endl;
        //   true_count++;
        // }
        // else {
        //   cout << "false" << endl;
        //   false_count++;
        // }

      }
    }
  }
  // cout << "true count" << true_count << endl;
  // cout << "false count" << false_count << endl;
}
  
