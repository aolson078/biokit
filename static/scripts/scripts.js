
function showSubRadios(reportId) {
    var subRadios = document.querySelectorAll('.sub-radios');
    subRadios.forEach(function (subRadio) {
        subRadio.style.display = 'none';
    });
    var subRadioId = 'sub-' + reportId;
    var subRadioElement = document.getElementById(subRadioId);
    if (subRadioElement) {
        subRadioElement.style.display = 'block';
    }
}




// --------------------------------------------------------------------------------------------------------------------
// These are for the drop down menus in the admin page to change username and password


document.getElementById('user-dropdown').addEventListener('change', function() {
    const selectedUserId = document.getElementById('user-dropdown').value;
    console.log('Selected User ID:', selectedUserId);});


const changeUsernameBtn = document.getElementById('changeUsernameBtn');
const usernameInput = document.getElementById('usernameInput');



usernameInput.addEventListener('change', (event) => {
    const newUsername = event.target.value;
    console.log(`New username: ${newUsername}`);
    });

document.getElementById('user-dropdown').addEventListener('change', function() {
    const selectedUserId = document.getElementById('user-dropdown').value;
    console.log('Selected User ID:', selectedUserId);});


function changeUsername() {
    const userDropdown = document.getElementById('user-dropdown');
    const user_id = userDropdown.value;
    const new_username = document.getElementById('usernameInput').value

    // AJAX request to Flask route
    fetch(`/change_username/${user_id}`, {
        method: 'PUT',
        headers: {
            'Content-Type': 'application/json',
        },
        body: JSON.stringify({ new_username })
    })
    .then(response => {
        if (response.ok) {
            console.log('Username changed successfully');
        } else {
            console.error('Error changing username')}
    })
}

function changePassword() {
    const userDropdown = document.getElementById('user-dropdown');
    const user_id = userDropdown.value;
    const new_password = document.getElementById('passwordInput').value
    console.log(new_password)
    // AJAX request to Flask route
    fetch(`/change_password/${user_id}`, {     // Get change password route
        method: 'PUT',
        headers: {
            'Content-Type': 'application/json',
        },
        body: JSON.stringify({ new_password })
    })
    .then(response => {
        if (response.ok) {
            console.log('Password changed successfully');
        } else {
            console.error('Error changing password')}
    })
}


// --------------------------------------------------------------------------------------------------------------------
// These are for the employee screen, to populate the cards

window.onload = function() {
//    localStorage.clear();
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


// Figure out if these are still needed!


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

// --------------------------------------------------------------------------------------------------------------------