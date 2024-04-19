//
//
//
//
//
//function showSubRadios(reportId) {
//    var subRadios = document.querySelectorAll('.sub-radios');
//    subRadios.forEach(function (subRadio) {
//        subRadio.style.display = 'none';
//    });
//    var subRadioId = 'sub-' + reportId;
//    var subRadioElement = document.getElementById(subRadioId);
//    if (subRadioElement) {
//        subRadioElement.style.display = 'block';
//    }
//}
//
//function compileReport() {
//        window.location.href = '/compile_report';
//    }
//
//function displayReport(reportId) {
//        window.location.href = '/display_report/' + reportId;
//    }
//
//// --------------------------------------------------------------------------------------------------------------------
//// These are for the drop down menus in the admin page and to change/create username and password
//
//// Puts user ID into Selected user on page
//window.onload = function() {
//    const params = new URLSearchParams(window.location.search);
//    const userId = params.get('user'); // get the user parameter from the URL
//
//    if (userId) {
//        const userIdSpan = document.querySelector('.user-id'); // select the span with the class "user-id"
//        if (userIdSpan) {
//            userIdSpan.textContent = userId; // set the text content of the span to the user ID
//        } else {
//            console.error('Could not find user ID span');
//        }
//    } else {
//        console.error('No user ID in URL');
//    }
//};
//
//
//
//document.getElementById('user-dropdown').addEventListener('change', function() {
//    const selectedUserId = document.getElementById('user-dropdown').value;
//    document.querySelector('.user-id').textContent = selectedUserId;
//    console.log('Selected User ID:', selectedUserId);});
//
//
//const changeUsernameBtn = document.getElementById('changeUsernameBtn');
//const usernameInput = document.getElementById('usernameInput');
//
//
//
//usernameInput.addEventListener('change', (event) => {
//    const newUsername = event.target.value;
//    console.log(`New username: ${newUsername}`);
//    });
//
//document.getElementById('user-dropdown').addEventListener('change', function() {
//    const selectedUserId = document.getElementById('user-dropdown').value;
//    console.log('Selected User ID:', selectedUserId);});
//
//// make delete user function (Will)
//function deleteUser() {
//
//}
//
//function changeUsername() {
//    const userDropdown = document.getElementById('user-dropdown');
//    const user_id = userDropdown.value;
//    const new_username = document.getElementById('usernameInput').value
//
//    // AJAX request to Flask route
//    fetch(`/change_username/${user_id}`, {
//        method: 'PUT',
//        headers: {
//            'Content-Type': 'application/json',
//        },
//        body: JSON.stringify({ new_username })
//    })
//    .then(response => {
//        if (response.ok) {
//            console.log('Username changed successfully');
//        } else {
//            console.error('Error changing username')}
//    })
//}
//
//function changePassword() {
//    const userDropdown = document.getElementById('user-dropdown');
//    const user_id = userDropdown.value;
//    const new_password = document.getElementById('passwordInput').value
//    console.log(new_password)
//    // AJAX request to Flask route
//    fetch(`/change_password/${user_id}`, {     // Get change password route
//        method: 'PUT',
//        headers: {
//            'Content-Type': 'application/json',
//        },
//        body: JSON.stringify({ new_password })
//    })
//    .then(response => {
//        if (response.ok) {
//            console.log('Password changed successfully');
//        } else {
//            console.error('Error changing password')}
//    })
//}
//
//
//// --------------------------------------------------------------------------------------------------------------------
//// These are for the employee screen, to populate the cards (Maybe take this out and display
//// cards on seperate route
//
//window.onload = function() {
////    localStorage.clear();
//    retrieveInputs();
//};
//
//
//function generateCard(input) {
//    let card = document.createElement("div");
//    card.innerHTML =
//        `
//        <div class="card">
//        <h3>Result:</h3>
//        // whats being put on the card
//        <p>${input}</p>
//        <a href = "http://supercoolwebsite.com">Click here!</a>
//        </div>
//        `;
//    return card;
//}
//
//
//function appendCard(card) {
//    let contentDiv = document.querySelector('.card-container');
//    contentDiv.appendChild(card);
//}
//
//
//function storeInput(input) {
//    let key = Date.now().toString();
//    localStorage.setItem(key, input);
//}
//
//
//function retrieveInputs() {
//    let keys = Object.keys(sessionStorage);
//    for (let key of keys) {
//        let input = sessionStorage.getItem(key);
//        if (input.trim() !== '') { // Check if input is not empty
//            let card = generateCard(input);
//            appendCard(card);
//        }
//    }
//}
