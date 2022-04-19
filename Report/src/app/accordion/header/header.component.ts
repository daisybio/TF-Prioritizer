import {Component, Input, OnInit} from '@angular/core';

@Component({
  selector: 'accordion-header',
  templateUrl: './header.component.html',
  styleUrls: ['./header.component.css']
})
export class HeaderComponent implements OnInit {

  @Input()
  title: string = "undefined";

  @Input()
  disabled: boolean = false;

  @Input()
  sub: boolean = false;

  constructor() {
  }

  ngOnInit(): void {
  }

  getClasses() {
    return {
      "light": this.sub
    }
  }
}
