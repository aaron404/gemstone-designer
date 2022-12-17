#[cfg(test)]

use proc_macro::TokenStream;
use quote::quote;

#[proc_macro]
pub fn get_uniform_location(tokens: proc_macro::TokenStream) -> proc_macro::TokenStream {
    //proc_macro::TokenStream::new()
    proc_macro::TokenStream::from(quote!{asdf})
}