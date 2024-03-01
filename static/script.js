//window.onload = function() {
//    localStorage.clear();
////    sessionStorage.clear();
//    var url = new URL(window.location.href);
//    var param = url.searchParams.get("parameter");
//    if (param) {
//        let decodedParam = decodeURIComponent(param);
//        let card = generateCard(decodedParam);
//        appendCard(card);
//        storeInput(decodeURIComponent(param));
//    }
//    retrieveInputs();
//};

window.onload = function() {
    localStorage.clear();
    retrieveInputs();
};


function generateCard(input) {
    let card = document.createElement("div");
    card.innerHTML =
        `
        <div class="card">
        <h3>Result:</h3>
        <p>${input}</p>
        <a href = "http://supercoolwebsite.com">Click here!</a>
        </div>
        `;
    return card;
}


function appendCard(card) {
    let contentDiv = document.querySelector('.card-container');
    contentDiv.appendChild(card);
}


function storeInput(input) {
    let key = Date.now().toString();
    localStorage.setItem(key, input);
}


function retrieveInputs() {
    let keys = Object.keys(sessionStorage);
    for (let key of keys) {
        let input = sessionStorage.getItem(key);
        if (input.trim() !== '') { // Check if input is not empty
            let card = generateCard(input);
            appendCard(card);
        }
    }
}
//function retrieveInputs() {
//    let keys = Object.keys(localStorage);
//    for (let key of keys) {
//        let input = localStorage.getItem(key);
//        if (input.trim() !== '') { // Check if input is not empty
//            let card = generateCard(input);
//            appendCard(card);
//        }
//    }
//}
//
//document.addEventListener("DOMContentLoaded", function() {
//  let form = document.querySelector("form[action='/search']");
//
//  form.addEventListener("submit", function(event) {
//    let input = document.querySelector("input[name='search']").value;
//
//      console.log(input);
//      storeInput(input);
//      let card = generateCard(input);
//      appendCard(card);
//      console.log("appended card:", card);
//      form.submit(); // This line submits the form programmatically
//
//  });
//});
