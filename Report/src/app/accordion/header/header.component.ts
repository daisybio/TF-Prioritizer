import {Component, Input, OnInit} from '@angular/core';

@Component({
  selector: 'accordion-header',
  templateUrl: './header.component.html',
  styleUrls: ['./header.component.css']
})
export class HeaderComponent implements OnInit {

  @Input()
    // @ts-ignore
  title: string;

  @Input()
    // @ts-ignore
  sub: boolean;

  constructor() {
  }

  ngOnInit(): void {
  }

  getClasses() {
    return {
      "sub": this.sub
    }
  }
}
